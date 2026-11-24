
"""
Command-line launcher for the POLLUX / HWO SNR calculator.
 
After `pip install -e .` (or `make install`) this is available as:
 
    pollux_etc                 # start the app and open http://127.0.0.1:8000
    pollux_etc --port 8080     # choose a port
    pollux_etc --no-browser    # do not open a browser tab
"""
import argparse
import os
import socket
import sys
from pathlib import Path
 
# Folder containing app3.py and the local modules (Telescope.py, ...).
# The original scripts use paths relative to src/, so we run from here.
APP_DIR = Path(__file__).resolve().parent
 
 
def _port_is_free(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.3)
        return s.connect_ex(("127.0.0.1", port)) != 0
 
 
def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="pollux_etc",
        description="Launch the POLLUX exposure-time/SNR calculator on localhost.",
    )
    parser.add_argument("--host", default="127.0.0.1",
                        help="interface to bind (default: 127.0.0.1)")
    parser.add_argument("-p", "--port", type=int, default=8000,
                        help="port to use (default: 8000; next free port is used if busy)")
    parser.add_argument("--no-browser", action="store_true",
                        help="do not open the browser automatically")
    args = parser.parse_args(argv)
 
    os.chdir(APP_DIR)
    sys.path.insert(0, str(APP_DIR))
 
    if not os.environ.get("PYSYN_CDBS"):
        print("WARNING: PYSYN_CDBS is not set; stsynphot reference data may not be found.",
              file=sys.stderr)
 
    port = args.port
    for _ in range(20):
        if _port_is_free(port):
            break
        print(f"Port {port} is busy, trying {port + 1}...")
        port += 1
    else:
        sys.exit(f"No free port found in {args.port}-{port}.")
 
    try:
        import app3 as app
    except ImportError as e:
        sys.exit(f"Could not start the app ({e}).\n"
                 "Install the dependencies with: pip install nicegui plotly numpy")
 
    print(f"POLLUX ETC running at http://{args.host}:{port}  (Ctrl+C to stop)")
    app.main(host=args.host, port=port, show=not args.no_browser, reload=False)
 
 
if __name__ == "__main__":
    main()
 

