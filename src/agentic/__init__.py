"""Post-ranking agent layer for CellLineSelector.

The modules in this package may explain and verify a deterministic result, but
they never calculate recommendation scores or change candidate order.
"""

from .output_agent import OutputAgent

__all__ = ["OutputAgent"]
