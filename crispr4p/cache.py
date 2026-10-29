"""Result cache."""

import logging
import os
import pickle
import tempfile


def ensure_cache_dir(cache_directory):
    """Create the cache directory if needed."""
    os.makedirs(cache_directory, exist_ok=True)


def cache_path(
    cache_directory,
    version,
    n_mismatch,
    chromosome,
    start,
    end,
    systematic_name=None,
):
    """Return a gene or coordinate cache path."""
    if systematic_name:
        basename = "%s_v%s_n%s.pickle" % (
            systematic_name,
            version,
            n_mismatch,
        )
    else:
        basename = "%s_%s_%s_v%s_n%s.pickle" % (
            chromosome,
            start,
            end,
            version,
            n_mismatch,
        )
    return os.path.join(cache_directory, basename)


def cache_exists(cache_path):
    return os.path.isfile(cache_path)


def load_cache(cache_path):
    """Load a cached design result."""
    with open(cache_path, "rb") as cache_file:
        payload = pickle.load(cache_file)
    return _normalize_result(payload)


def save_cache(cache_path, result):
    """Publish a complete pickle without exposing a partial write to readers."""
    temporary_path = None
    try:
        data = pickle.dumps(_normalize_result(result), protocol=-1)
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=os.path.dirname(cache_path) or ".",
            prefix=".cache-", suffix=".tmp", delete=False,
        ) as cache_file:
            temporary_path = cache_file.name
            cache_file.write(data)
        os.replace(temporary_path, cache_path)
    except Exception as error:
        logging.getLogger(__name__).warning(
            "Could not save cache %s: %s", cache_path, error
        )
    finally:
        if temporary_path is not None:
            try:
                os.unlink(temporary_path)
            except FileNotFoundError:
                pass
            except OSError as error:
                logging.getLogger(__name__).warning(
                    "Could not remove cache temporary file %s: %s", temporary_path, error
                )


def drop_cache(cache_path):
    """Remove a cached result."""
    if os.path.isfile(cache_path):
        os.remove(cache_path)


def get_cached(cache_path, compute, is_cached=cache_exists):
    """Load a result or compute it."""
    if is_cached(cache_path):
        try:
            return load_cache(cache_path)
        except Exception:
            # A different worker may already be replacing this file; don't unlink it.
            pass
    result = _normalize_result(compute())
    save_cache(cache_path, result)
    return result


def _normalize_result(result):
    table, hr_dna, checking_primers, guide_matches = result
    return table, hr_dna, checking_primers, guide_matches
