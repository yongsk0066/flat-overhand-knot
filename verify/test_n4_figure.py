"""Check visibility at the coincident hexagon sides and the four tail arrows in Figure 9."""

from pathlib import Path
import unittest
import xml.etree.ElementTree as ET

from shapely.geometry import LineString
from shapely.ops import unary_union


def points(element):
    return [tuple(map(float, p.split(","))) for p in element.get("points").split()]


class TightFigure(unittest.TestCase):
    def setUp(self):
        path = Path(__file__).resolve().parents[1] / "paper/figs/fig41-tight.svg"
        self.svg = ET.parse(path).getroot()
        self.views = [points(p) for p in self.svg.findall("{*}polygon")
                      if p.get("stroke") == "#0a0a0a"]
        self.assertEqual([len(v) for v in self.views], [6, 6])

    def test_tail_edges_under_the_core(self):
        corners = self.views[1]
        for color, a, b in [("#c1121f", 1, 2), ("#8ecae6", 4, 5)]:
            side = LineString([corners[a], corners[b]])
            hidden, visible = [], []
            for line in self.svg.findall("{*}polyline"):
                if line.get("stroke") == color:
                    group = hidden if "stroke-dasharray" in line.attrib else visible
                    group.append(LineString(points(line)))
            self.assertEqual(side.difference(unary_union(hidden)).length, 0, color)
            self.assertEqual(side.intersection(unary_union(visible)).length, 0, color)

    def test_tail_arrows_in_both_views(self):
        width, height = float(self.svg.get("width")), float(self.svg.get("height"))
        arrows = [tuple(float(line.get(k)) for k in ("x1", "y1", "x2", "y2"))
                  for line in self.svg.findall("{*}line") if "marker-end" in line.attrib]
        self.assertEqual(len(arrows), 4)
        shift = self.views[1][0][0] - self.views[0][0][0]
        for i, (x1, y1, x2, y2) in enumerate(arrows):
            left, right = (0, width / 2) if i < 2 else (width / 2, width)
            self.assertTrue(left < min(x1, x2) <= max(x1, x2) < right and 0 < y1 == y2 < height)
            self.assertEqual(x2 < x1, i % 2 == 0)
            if i < 2:
                self.assertEqual(arrows[i + 2], (x1 + shift, y1, x2 + shift, y2))


if __name__ == "__main__":
    unittest.main()
