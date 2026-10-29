#!/usr/bin/env python3
# Web entry points and local preview server.

import json
import os
import urllib.parse
import sys
from functools import lru_cache
from html import escape
from http import HTTPStatus
from io import BytesIO

from crispr4p.service import (
    AmbiguousGeneNameError,
    Crispr4pService,
    GeneNameNotFoundError,
    OligoLengthError,
    PrimerNotFoundError,
)
from crispr4p.web_views import (
    render_design,
    render_ambiguous_gene_error,
    render_error,
    render_gene_error,
    render_oligo,
    render_oligo_error,
    render_query_error,
)


PORT = 8080
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
MAX_FORM_SIZE = 64 * 1024
DONOR_ARM_LENGTH = 80
CASSETTE_LENGTH = 23
PRIMER_WINDOW = 300
NEB_LINKS = {
    "PacI": (
        "https://www.neb.com/en-us/products/r0547-paci",
        "https://nebcloner.neb.com/#!/protocol/re/single/PacI",
    ),
    "SwaI": (
        "https://www.neb.com/en-us/products/r0604-swai",
        "https://nebcloner.neb.com/#!/protocol/re/single/SwaI",
    ),
}


def primer_payload(pair, checks=None):
    if pair is None:
        return None
    payload = {
        "forward": pair.forward,
        "reverse": pair.reverse,
        "forward_tm": pair.forward_tm,
        "reverse_tm": pair.reverse_tm,
        "wt_product_size": pair.wt_product_size,
        "disrupted_product_size": pair.disrupted_product_size,
    }
    if checks is not None:
        payload["left_junction"] = {
            "forward": checks.left.forward,
            "reverse": checks.left.reverse,
            "forward_tm": checks.left.forward_tm,
            "reverse_tm": checks.left.reverse_tm,
            "product_size": checks.left.product_size,
        }
        payload["right_junction"] = {
            "forward": checks.right.forward,
            "reverse": checks.right.reverse,
            "forward_tm": checks.right.forward_tm,
            "reverse_tm": checks.right.reverse_tm,
            "product_size": checks.right.product_size,
        }
    return payload


@lru_cache(maxsize=1)
def create_service():
    """Create the shared application service."""
    return Crispr4pService.from_project_data(
        precomputed_folder=os.path.join(PROJECT_DIR, "precomputed"),
    )


class WebHandler:
    def do_GET(self):
        parsed_path = urllib.parse.urlparse(self.path)
        path = parsed_path.path
        query = urllib.parse.parse_qs(parsed_path.query)

        # Query endpoints also work when Apache mounts only the CGI script.
        if path in ("/", "/index.html", "/webapp.py"):
            path = "/" + query.get("endpoint", [""])[0]

        if path == "/":
            self.serve_form()
        elif path == "/insertion-primers":
            self.serve_insertion_primers(query)
        elif path == "/cassette-options":
            self.serve_cassette_options(query)
        elif path == "/css":
            self.serve_css("crispr4p.css")
        elif path.startswith("/css/"):
            self.serve_css(path)
        else:
            self.send_error(404, "File not found")

    def do_POST(self):
        parsed_path = urllib.parse.urlparse(self.path)
        path = parsed_path.path

        if path == "/webapp.py" or path == "/":
            self.process_post()
        else:
            self.send_error(404, "File not found")

    def serve_form(self, result_content=""):
        try:
            template_path = os.path.join(PROJECT_DIR, 'template/bahler_template.html')
            with open(template_path, 'r', encoding='utf-8') as fh:
                template_file = fh.read()
            
            # Insert the result into the page template.
            rendered = template_file % result_content
            
            self.send_response(200)
            self.send_header("Content-type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(rendered.encode('utf-8'))
        except IOError as err:
            self.send_error(500, f"Error reading templates: {err}")

    def serve_css(self, path):
        filename = os.path.basename(path)
        css_file = os.path.join(PROJECT_DIR, "css", filename)
        if os.path.exists(css_file):
            try:
                with open(css_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                self.send_response(200)
                self.send_header("Content-type", "text/css")
                self.end_headers()
                self.wfile.write(content.encode('utf-8'))
            except IOError:
                self.send_error(500, "Error reading CSS file")
        else:
            self.send_error(404, "CSS File not found")

    def serve_json(self, payload, status=200):
        content = json.dumps(payload, separators=(",", ":")).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-type", "application/json; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(content)))
        self.end_headers()
        self.wfile.write(content)

    def serve_insertion_primers(self, query):
        try:
            chromosome = query["chromosome"][0].strip()
            cut = (
                int(query["cut_left"][0]),
                int(query["cut_right"][0]),
            )
            cassette_id = query.get("cassette_id", [None])[0]
            coding_strand = query.get("coding_strand", [None])[0]
            if (cassette_id is None) != (coding_strand is None):
                raise ValueError(
                    "cassette_id and coding_strand must be provided together"
                )

            service = create_service()
            checks = None
            if cassette_id is None:
                pair = service.insertion_primers(
                    chromosome,
                    cut,
                    arm_length=DONOR_ARM_LENGTH,
                    insert_length=CASSETTE_LENGTH,
                    window=PRIMER_WINDOW,
                )
            else:
                checks = service.insertion_checks(
                    chromosome,
                    cut,
                    int(cassette_id),
                    coding_strand,
                    arm_length=DONOR_ARM_LENGTH,
                    window=PRIMER_WINDOW,
                )
                pair = checks.spanning
        except PrimerNotFoundError:
            self.serve_json(
                {"error": "No insertion-checking primer pair was found."},
                status=422,
            )
            return
        except (KeyError, IndexError, ValueError) as error:
            self.serve_json({"error": str(error)}, status=400)
            return

        self.serve_json(primer_payload(pair, checks))

    def serve_cassette_options(self, query):
        try:
            chromosome = query["chromosome"][0].strip()
            cut = (
                int(query["cut_left"][0]),
                int(query["cut_right"][0]),
            )
            guide = query["guide"][0].strip()
            cassette_id = int(query["cassette_id"][0])
            coding_strand = query["coding_strand"][0]
            options = create_service().cassette_options(
                chromosome,
                cut,
                guide,
                cassette_id,
                coding_strand,
                arm_length=DONOR_ARM_LENGTH,
                window=PRIMER_WINDOW,
            )
        except (KeyError, IndexError, ValueError) as error:
            self.serve_json({"error": str(error)}, status=400)
            return

        rows = []
        for option in options:
            item = option.cassette_format
            digest = option.digest
            product_url, protocol_url = NEB_LINKS.get(
                item.enzyme,
                (None, None),
            )
            rows.append({
                "id": item.id,
                "label": item.label,
                "length": item.length,
                "enzyme": item.enzyme,
                "site": item.site,
                "product_url": product_url,
                "protocol_url": protocol_url,
                "available": option.available,
                "coding_sequence": option.coding_sequence,
                "insert": option.insert,
                "primers": primer_payload(option.spanning, option.checks),
                "digest": (
                    {
                        "wt_site_count": len(digest.wt_sites),
                        "edited_site_count": len(digest.edited_sites),
                        "fragments": digest.fragments,
                    }
                    if digest is not None
                    else None
                ),
            })
        self.serve_json({"formats": rows})

    def process_post(self):
        content_type = self.headers.get(
            "Content-Type", "application/x-www-form-urlencoded"
        ).split(";", 1)[0].strip().lower()
        if content_type != "application/x-www-form-urlencoded":
            self.send_error(415, "Expected a URL-encoded form")
            return
        try:
            content_length = int(self.headers.get('Content-Length', 0))
            if content_length < 0:
                raise ValueError
        except ValueError:
            self.send_error(400, "Invalid Content-Length")
            return
        if content_length > MAX_FORM_SIZE:
            self.send_error(413, "Form is too large")
            return
        try:
            body = self.rfile.read(content_length)
            if len(body) != content_length:
                raise ValueError
            post_data = body.decode('utf-8')
        except ValueError:
            self.send_error(400, "Invalid form body")
            return
        params = urllib.parse.parse_qs(post_data)

        name = params.get('name', [None])[0]
        chromosome = params.get('chromosome', [None])[0]
        coor_lower = params.get('coor_lower', [None])[0]
        coor_upper = params.get('coor_upper', [None])[0]
        oligo_sequence = params.get('oligo_sequence', [None])[0]
        oligo_mismatch_str = params.get('oligo_mismatch', ['0'])[0]

        name = name.strip() if name else None
        chromosome = chromosome.strip() if chromosome else None
        coor_lower = coor_lower.strip() if coor_lower else None
        coor_upper = coor_upper.strip() if coor_upper else None
        oligo_sequence = oligo_sequence.strip().upper() if oligo_sequence else None

        result_html = ""

        try:
            if oligo_sequence:
                try:
                    mismatches = int(oligo_mismatch_str)
                except ValueError:
                    mismatches = 0
                result_html = self.run_oligo_model(oligo_sequence, mismatches)
            elif name or (chromosome and coor_lower and coor_upper):
                result_html = self.run_design_model(name, chromosome, coor_lower, coor_upper)
            else:
                result_html = render_query_error()
        except GeneNameNotFoundError as error:
            result_html = render_gene_error(error.query)
        except AmbiguousGeneNameError as error:
            result_html = render_ambiguous_gene_error(error)
        except Exception as e:
            result_html = render_error(e)

        self.serve_form(result_html)

    def run_oligo_model(self, oligo_seq, mismatches):
        try:
            result = create_service().analyze_oligo(
                oligo_seq,
                n_mismatch=mismatches,
            )
        except OligoLengthError as error:
            return render_oligo_error(error.sequence_length)

        return render_oligo(result)

    def run_design_model(self, name, chromosome, coor_lower, coor_upper):
        service = create_service()
        if name is not None:
            result = service.design_gene(name, n_mismatch=0)
        else:
            result = service.design_region(
                chromosome,
                coor_lower,
                coor_upper,
                n_mismatch=0,
            )
        guide_annotations = service.annotate_guides(result.guides)
        cassette_choices = service.cassette_choices(
            result.guides,
            guide_annotations,
            result.name,
        )
        disruption_donors = service.disruption_donors(
            result.guides,
            guide_annotations,
            cassette_choices,
            DONOR_ARM_LENGTH,
            result.name,
        )
        restoration_donors = service.restoration_donors(
            result.guides,
            guide_annotations,
            DONOR_ARM_LENGTH,
            result.name,
        )

        with open(os.path.join(PROJECT_DIR, 'template/container_table.html'), 'r', encoding='utf-8') as fh:
            template_file = fh.read()

        return render_design(
            result,
            guide_annotations,
            template_file,
            cassette_choices=cassette_choices,
            disruption_donors=disruption_donors,
            restoration_donors=restoration_donors,
        )


if "GATEWAY_INTERFACE" not in os.environ:
    # Apache handles HTTP itself; this is only needed for local preview.
    from http.server import HTTPServer, BaseHTTPRequestHandler

    class CRISPR4PHandler(WebHandler, BaseHTTPRequestHandler):
        pass


class RequestHandler(WebHandler):
    def __init__(self, environ, stdin, stdout):
        self.command = environ.get("REQUEST_METHOD", "GET")
        self.path = environ.get("PATH_INFO") or "/"
        self.path += "?" + environ.get("QUERY_STRING", "")
        self.headers = {
            "Content-Length": environ.get("CONTENT_LENGTH", "0"),
            "Content-Type": environ.get(
                "CONTENT_TYPE", "application/x-www-form-urlencoded"
            ),
        }
        self.rfile = stdin
        self.wfile = stdout

    def run(self):
        if self.command == "GET":
            self.do_GET()
        elif self.command == "POST":
            self.do_POST()
        else:
            self.send_error(405, "Method not allowed")

    def send_error(self, status, message):
        content = f"<h1>{status}</h1><p>{escape(message)}</p>".encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(content)))
        if status == 405:
            self.send_header("Allow", "GET, POST")
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(content)


class CGIHandler(RequestHandler):
    def __init__(self, environ, stdin, stdout):
        super().__init__(environ, stdin, stdout)
        self.path = "/webapp.py" + environ.get("PATH_INFO", "")
        self.path += "?" + environ.get("QUERY_STRING", "")

    def send_response(self, status):
        self.send_header("Status", f"{status} {HTTPStatus(status).phrase}")

    def send_header(self, name, value):
        self.wfile.write(f"{name}: {value}\r\n".encode("latin-1"))

    def end_headers(self):
        self.wfile.write(b"\r\n")


class WSGIHandler(RequestHandler):
    def __init__(self, environ):
        super().__init__(environ, environ["wsgi.input"], BytesIO())
        self.response_headers = []

    def send_response(self, status):
        self.status = f"{status} {HTTPStatus(status).phrase}"

    def send_header(self, name, value):
        self.response_headers.append((name, str(value)))

    def end_headers(self):
        pass


def application(environ, start_response):
    """WSGI entry point; Gunicorn workers retain their application service."""
    handler = WSGIHandler(environ)
    handler.run()
    body = handler.wfile.getvalue()
    if not any(name.lower() == "content-length" for name, _ in handler.response_headers):
        handler.send_header("Content-Length", len(body))
    start_response(handler.status, handler.response_headers)
    return [body]


def main():
    if "GATEWAY_INTERFACE" in os.environ:
        CGIHandler(os.environ, sys.stdin.buffer, sys.stdout.buffer).run()
        return

    server_address = ('', PORT)
    httpd = HTTPServer(server_address, CRISPR4PHandler)
    print(f"Starting local server on http://localhost:{PORT} ...")
    print("Press Ctrl+C to stop.")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server.")
        sys.exit(0)


if __name__ == "__main__":
    main()
