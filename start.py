"""Single entry point: launch the UI by default, or use ``--cli``."""

from __future__ import annotations

import sys


def main() -> None:
    if "--cli" in sys.argv:
        sys.argv.remove("--cli")
        from dynamic_cellline_selector_gene_protein import main as cli_main

        cli_main()
        return

    from api_server import run_server

    run_server()


if __name__ == "__main__":
    main()
