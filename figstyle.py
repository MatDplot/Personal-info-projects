"""Optional house style for the figures.

If a private `_style/` folder sits next to this file (it is not published), the
figures get the look of my thesis: TeX Gyre Heros and Computer Modern, text set
by lualatex through the pgf backend, numbers and units by siunitx, legends without
a frame, their text outlined in white (contour package) where it crosses a curve.
Without it, nothing is required beyond matplotlib: labels fall back to mathtext.

    import figstyle
    figstyle.apply()
    ax.set_xlabel(f"$E$ {figstyle.unit(r'\\mega\\electronvolt', 'MeV')}")
    ax.set_title(f"$L = {figstyle.num(48)}$")
"""
import sys
from pathlib import Path

_PRIVATE = Path(__file__).resolve().parent / "_style"
TEX = False
STYLE = None


def apply():
    """Install the private style if it is there. Returns its module, or None."""
    global TEX, STYLE
    import matplotlib as mpl
    mpl.rcParams["figure.constrained_layout.use"] = True
    if not _PRIVATE.is_dir():
        return None
    sys.path.insert(0, str(_PRIVATE))
    try:
        import plotstyle
        import thesis_tex
    except ImportError:
        return None
    TEX = thesis_tex.enable(quiet=True)
    plotstyle.apply()
    if TEX:
        _fix_minus()
        mpl.rcParams["axes.formatter.use_mathtext"] = False   # tick labels in the text font
    mpl.rcParams.update({"legend.frameon": False, "figure.constrained_layout.use": True})
    STYLE = plotstyle
    _thesis_geometry_and_ticks()
    _contoured_legends()
    return plotstyle


def unit(si, plain):
    """A unit in brackets, the brackets inside the maths as in the thesis: $[\\unit{...}]$.

    siunitx macros when LaTeX renders the text, mathtext otherwise.

    unit(r"\\giga\\electronvolt", r"GeV"), unit(r"\\metre\\per\\second", r"m\\,s^{-1}")
    """
    return rf"$[\unit{{{si}}}]$" if TEX else rf"$[\mathrm{{{plain}}}]$"


def num(value):
    """A number set by siunitx, \\num{value}, under LaTeX (in text or maths); the plain number otherwise."""
    return rf"\num{{{value}}}" if TEX else f"{value}"


def halo(text):
    """A free text label lifted off the curves behind it: \\contour{white}{...} under LaTeX."""
    if not TEX:
        return text
    import thesis_text
    return thesis_text.contoured(text)


def cmap(name, fallback):
    """A colour map of the private style (CMAP_DIVERGING, CMAP_BLUE, ...) or a matplotlib one."""
    return getattr(STYLE, name, fallback)


def qty(value, si, plain):
    """A value with its unit: \\qty{value}{si} under LaTeX, mathtext otherwise."""
    return rf"\qty{{{value}}}{{{si}}}" if TEX else rf"${value}\,\mathrm{{{plain}}}$"


def _fix_minus():
    """U+2212 as \\textminus in text but a plain - inside $...$, where \\textminus prints nothing."""
    import re
    from matplotlib.backends import backend_pgf as pgf
    inner = pgf._tex_escape

    def escape(text):
        parts = re.split(r"(?<!\\)\$", text)
        for i in range(1, len(parts), 2):
            parts[i] = parts[i].replace("\N{MINUS SIGN}", "-")
        return inner("$".join(parts))

    pgf._tex_escape = escape


def _contoured_legends():
    """Every legend gets the thesis halo: \\contour behind its text under pgf, a path-effect halo on Agg (GIF frames)."""
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    import thesis_text

    def wrap(make):
        def legend(self, *args, **kwargs):
            leg = make(self, *args, **kwargs)
            fig = self if isinstance(self, Figure) else self.figure
            if type(fig.canvas).__name__ == "FigureCanvasPgf":
                # siunitx does not survive \contour, and a path-effect halo would send the label to
                # mathtext: those labels stay without a halo
                for text in leg.get_texts():
                    text.set_text(thesis_text.contoured(text.get_text()))
            else:
                thesis_text.halo(leg.get_texts())
            return leg
        return legend

    Axes.legend = wrap(Axes.legend)
    Figure.legend = wrap(Figure.legend)


def _thesis_geometry_and_ticks():
    """At save time: thesis text width (16 cm) and the thesis tick formatter, as in Thesis-plot."""
    import matplotlib.ticker as ticker
    from matplotlib.figure import Figure
    from thesis_ticks import format_tick

    width = STYLE.THESIS_TEXT_WIDTH_PT / 72.27
    plain = (ticker.ScalarFormatter, ticker.LogFormatterSciNotation)
    save = Figure.savefig

    def savefig(self, *args, **kwargs):
        if type(self.canvas).__name__ == "FigureCanvasPgf":   # not the Agg-drawn GIF frames
            w, h = self.get_size_inches()
            self.set_size_inches(width, h * width / w)
            for ax in self.axes:
                for axis in (ax.xaxis, ax.yaxis):
                    if isinstance(axis.get_major_formatter(), plain):
                        axis.set_major_formatter(ticker.FuncFormatter(lambda v, _p, a=axis: format_tick(v, a)))
        return save(self, *args, **kwargs)

    Figure.savefig = savefig
