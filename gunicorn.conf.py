bind = "127.0.0.1:8080"
workers = 4
# Keep Primer3 calls and mutable request state within one request per process.
worker_class = "sync"
threads = 1
timeout = 120
accesslog = "-"
errorlog = "-"
