"""
pyPLP_simulator/web_server.py: Lightweight local web server for pyPLP animated reaction simulator.
"""

import os
import sys
import webbrowser
import socket
from http.server import SimpleHTTPRequestHandler, HTTPServer
import threading
import time

def find_free_port(start_port=8080):
    port = start_port
    while port < start_port + 100:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(("127.0.0.1", port))
                return port
        except OSError:
            port += 1
    return start_port

def run_web_simulator(port=None, open_browser=True):
    web_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")
    html_file = os.path.join(web_dir, "web_simulator.html")
    
    if not os.path.exists(html_file):
        print(f"Error: Could not find web_simulator.html at {html_file}", file=sys.stderr)
        sys.exit(1)

    if port is None:
        port = find_free_port(8080)

    class CustomHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=web_dir, **kwargs)

        def do_GET(self):
            if self.path == "/" or self.path == "":
                self.path = "/web_simulator.html"
            return super().do_GET()

        def log_message(self, format, *args):
            # Suppress noisy standard request logging
            pass

    server = HTTPServer(("127.0.0.1", port), CustomHandler)
    url = f"http://127.0.0.1:{port}/"

    print("=" * 80)
    print("pyPLP Interactive Animated Molecular Simulator Server")
    print("=" * 80)
    print(f"Server running at: {url}")
    print("Serving interactive canvas animation for PLP target hybridization,")
    print("ligation nick verification, phi29 RCA rolling circle, and fluorescence.")
    print("Press Ctrl+C to stop the server.")
    print("=" * 80)

    if open_browser:
        def open_tab():
            time.sleep(0.5)
            webbrowser.open(url)
        threading.Thread(target=open_tab, daemon=True).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping web simulator server.")
        server.server_close()

if __name__ == "__main__":
    port_arg = int(sys.argv[1]) if len(sys.argv) > 1 else None
    run_web_simulator(port=port_arg)
