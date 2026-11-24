"""
NiceGUI front end for the HWO  + POLLUX SNR calculator.

Normally started with the `pollux_etc` command. For development:

    cd src
    python app.py --reload

Requires: nicegui>=2.0 plotly numpy (plus the local modules get_pollux_spectra,
pollux_snr_util, Telescope, planet_sim).
"""
import argparse
 
import numpy as np
import plotly.graph_objects as go
from nicegui import run as ng_run
from nicegui import ui
from plotly.subplots import make_subplots
 
import get_pollux_spectra
import pollux_snr_util as psu
import Telescope as T
 
PARSEC2AU = 206264.80624538
 
# Colour palette
C_MID = "#8D7E9D"    # muted purple
C_BG = "#F9F8FA"     # near white
C_DARK = "#363D8D"   # deep indigo
C_LIGHT = "#BDBBCA"  # light lavender grey
C_BLUE = "#1E4C96"   # blue
 
# Instrument modes for each observing type (first entry in each list is not special;
# DEFAULT_MODE gives the one selected when the type is switched).
MODES = {
    "Spectroscopy": ["FUV","MUV", "NUV", "OPT", "NIR"],
    "Polarimetry": ["FUVPOL", "MUVPOL", "NUVPOL", "OPTPOL", "NIRPOL"],
}
DEFAULT_MODE = {"Spectroscopy": "MUV", "Polarimetry": "MUVPOL"}
# Sources come from get_pollux_spectra. "Templates" are observed/synthetic spectra used
# as they are; "Models" are spectral-type models whose flux is scaled to a distance.
# Names are passed to get_pollux_spectra.add_spectrum_to_library() unchanged.
TEMPLATES = [
    "Procyon", "Tau Ceti", "HD 120411", "V Eps Eri", "WASP-33", "TRAPPIST-1",
    "GJ 832", "GJ 729", "HD 40307", "Classical T Tauri",
    "G191B2B (WD)", "GD71 (WD)", "GD153 (WD)", "QSO", "10 Myr Starburst",
    "Galaxy with f_esc, HI=1, HeI=1", "Galaxy with f_esc, HI=0.001, HeI=1",
    "Orion Nebula", "Starburst, No Dust", "Starburst, E(B-V) = 0.6",
]
MODELS = ["M7", "M3", "M1", "K2", "K2V", "G8", "G2V", "G0V", "F5", "A8"]
SOURCES = {"Model": MODELS, "Template": TEMPLATES}
DEFAULT_SOURCE = {"Model": "K2", "Template": "Procyon"}
 
 
# ----------------------------------------------------------------------------
# Core calculation (the original script, wrapped in a function)
# ----------------------------------------------------------------------------
def run_snr(aperture, temperature, diff_limit, mode, extime, template, is_model, distance):
    luvoir = T.Telescope(aperture, temperature, diff_limit)
    pollux = T.Spectropolarimeter()
    pollux.set_mode(mode)
 
    spec_dict = get_pollux_spectra.add_spectrum_to_library(template)
    if template not in spec_dict:
        raise KeyError(
            f"Source '{template}' was not found in the spectral library. "
            f"Available keys: {', '.join(map(str, spec_dict))}"
        )
    spec = spec_dict[template]
    wave = np.asarray(spec.wave)
    flux = np.array(spec.flux, dtype=float)  # copy: never mutate the library
 
    used_distance = None
    if is_model:
        flux = flux / (distance * PARSEC2AU**2)  # scale model flux to Earth
        spec.flux = flux
        spec.distance = distance
        used_distance = distance
 
    snr, source_counts, bef_int, b_counts = psu.simulate_exposure(
        luvoir, pollux, wave, flux, extime
    )
    snr = np.asarray(snr)
    source_counts = np.asarray(source_counts)
    bef_int = np.asarray(bef_int)
 
    st = psu.find_nearest(wave, pollux.lambda_range[0])
    en = psu.find_nearest(wave, pollux.lambda_range[1])
 
    warning = None
    if wave[-1] < pollux.lambda_range[1]:
        warning = (
            f"For mode {mode}, the source spectrum does not cover the "
            "entire instrument bandwidth."
        )
 
    header = psu.text_head_generator(
        pollux, template, **({"distance": used_distance} if used_distance is not None else {})
    )
 
    return dict(
        pollux=pollux, template=template, mode=mode, wave=wave, flux=flux,
        snr=snr, source_counts=source_counts, bef_int=bef_int, st=st, en=en,
        header=header, warning=warning, aperture=aperture,
    )
 
 
def empty_figure(message="Press 'Calculate S/N Ratio' to see results"):
    fig = go.Figure()
    fig.update_layout(
        template="plotly_white", height=720, paper_bgcolor=C_BG, plot_bgcolor=C_BG,
        font=dict(color=C_DARK),
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        annotations=[dict(text=message, showarrow=False, font=dict(size=16))],
    )
    return fig
 
 
def make_figure(r):
    pollux, wave, st, en = r["pollux"], r["wave"], r["st"], r["en"]
    lo, hi = pollux.lambda_range
    aeff = pollux.aeff * np.pi * (r["aperture"] * 1.0e2 / 2.0) ** 2
 
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{}, {"rowspan": 2}], [{}, None]],
        column_widths=[0.62, 0.38],
        horizontal_spacing=0.09, vertical_spacing=0.14,
        subplot_titles=(
            f"Flux vs Wavelength ({r['template']})",
            "Effective area",
            f"S/N Ratio vs Wavelength ({r['template']}, {pollux.name} mode {pollux.mode_name})",
        ),
    )
    fig.add_trace(go.Scattergl(x=wave, y=r["flux"], mode="lines", name="Flux",
                               line=dict(color="darkred", width=1.5)), row=1, col=1)
    fig.add_trace(go.Scattergl(x=wave, y=r["snr"], mode="lines", name="S/N",
                               line=dict(color="darkgreen", width=1.5)), row=2, col=1)
    fig.add_trace(go.Scattergl(x=np.asarray(pollux.wave), y=np.asarray(aeff), mode="lines",
                               name="A_eff", line=dict(color=C_MID, width=1.5)), row=1, col=2)
 
    # Flux and SNR panels share the x axis, limited to the instrument band
    fig.update_xaxes(range=[lo, hi], row=1, col=1)
    fig.update_xaxes(range=[lo, hi], matches="x", row=2, col=1)
    fig.update_yaxes(range=[0, float(np.max(r["flux"][st:en])) * 1.05], row=1, col=1)
    fig.update_xaxes(title_text="Wavelength [Å]")
    fig.update_yaxes(title_text="Flux [erg cm⁻² s⁻¹ Å⁻¹]", exponentformat="e", row=1, col=1)
    fig.update_yaxes(title_text="S/N per resel", row=2, col=1)
    fig.update_yaxes(title_text="Effective area [cm²]", exponentformat="e", row=1, col=2)
    fig.update_layout(
        height=720, showlegend=False, template="plotly_white",
        hovermode="x unified", margin=dict(l=70, r=20, t=50, b=50),
        paper_bgcolor=C_BG, plot_bgcolor=C_BG, font=dict(color=C_DARK),
    )
    fig.update_xaxes(gridcolor=C_LIGHT, linecolor=C_MID, zeroline=False)
    fig.update_yaxes(gridcolor=C_LIGHT, linecolor=C_MID, zeroline=False)
    return fig
 
 
def snr_table_bytes(r):
    """Same text table the original script writes to output/."""
    import io
 
    st, en = r["st"], r["en"]
    buf = io.StringIO()
    np.savetxt(
        buf,
        np.c_[r["wave"][st:en], r["flux"][st:en], r["bef_int"][st:en],
              r["source_counts"][st:en], r["snr"][st:en]],
        header=r["header"],
        fmt=["% 9.3f", "% 12.5e", "% 12.5e", "% 12.5e", "% 12.5e"],
    )
    return buf.getvalue().encode()
 
 
# ----------------------------------------------------------------------------
# UI (one independent state per browser session)
# ----------------------------------------------------------------------------
@ui.page("/")
def index():
    state = {"result": None}
 
    ui.page_title("POLLUX ETC")
    ui.colors(primary=C_DARK, secondary=C_MID, accent=C_BLUE, positive=C_BLUE, info=C_MID)
    ui.query("body").style(f"background-color: {C_BG}")
 
    with ui.header().classes("items-center"):
        ui.label("HWO + POLLUX: S/N Ratio calculator").classes("text-h6")
 
    with ui.left_drawer(fixed=True, bordered=True).props("width=330").classes("gap-1") \
            .style("background-color: rgba(189, 187, 202, 0.35)"):
        ui.label("Instrument").classes("text-subtitle1 text-bold")
        obs_type = ui.toggle(list(MODES), value="Spectroscopy").props("spread no-caps").classes("w-full")
        mode = ui.select(MODES["Spectroscopy"], value=DEFAULT_MODE["Spectroscopy"],
                         label="POLLUX mode").classes("w-full")
 
        def on_type_change(e):
            mode.set_options(MODES[e.value], value=DEFAULT_MODE[e.value])
 
        obs_type.on_value_change(on_type_change)
        extime = ui.number("Exposure time [hours]", value=1.0, min=0.001, step=0.5,
                           format="%.3f").classes("w-full")
 
        ui.separator()
        ui.label("Source").classes("text-subtitle1 text-bold")
        source_type = ui.toggle(list(SOURCES), value="Model").props("spread no-caps").classes("w-full")
        template = ui.select(SOURCES["Model"], value=DEFAULT_SOURCE["Model"], label="Source",
                             with_input=True).classes("w-full")
        distance = ui.number("Distance [pc]", value=10, min=0.01).classes("w-full")
        # Models are scaled to a distance; templates are used as provided
        distance.bind_visibility_from(source_type, "value", backward=lambda v: v == "Model")
 
        def is_model():
            return source_type.value == "Model"
 
        def on_source_change(e):
            template.set_options(SOURCES[e.value], value=DEFAULT_SOURCE[e.value])
 
        source_type.on_value_change(on_source_change)
 
        ui.separator()
        with ui.expansion("Telescope (HWO)", icon="settings").classes("w-full"):
            aperture = ui.number("Aperture [m]", value=8.0, min=0.1).classes("w-full")
            temperature = ui.number("Temperature [K]", value=280.0, min=1).classes("w-full")
            diff_limit = ui.number("Diffraction limit [nm]", value=500.0, min=1).classes("w-full")
 
        run_btn = ui.button("Calculate S/N Ratio", icon="calculate").classes("w-full mt-2")
        dl_txt = ui.button("Download S/N Ratio table (.txt)", icon="download").props("outline") \
            .classes("w-full mt-2")
        dl_html = ui.button("Download interactive plot (.html)", icon="download").props("outline") \
            .classes("w-full mt-2")
        dl_txt.disable()
        dl_html.disable()
 
    with ui.column().classes("w-full"):
        status = ui.column().classes("w-full gap-1")
        with ui.card().classes("w-full"):
            plot = ui.plotly(empty_figure()).classes("w-full")
        with status:
            ui.label("Set the parameters and press 'Calculate S/N Ratio'.").classes(
                "rounded p-3 w-full").style(f"background-color: {C_LIGHT}; color: {C_DARK}")
 
    def show_message(text, kind):
        styles = {
            "ok": f"background-color: {C_LIGHT}; color: {C_DARK}",
            "warn": f"background-color: {C_MID}; color: {C_BG}",
            "error": "background-color: #F8D7DA; color: #7A1C25",  # errors keep a red tone
        }
        with status:
            ui.label(text).classes("rounded p-3 w-full").style(styles[kind])
 
    async def compute():
        values = [aperture.value, temperature.value, diff_limit.value, extime.value]
        if any(v is None for v in values) or (is_model() and distance.value is None):
            ui.notify("Please fill in all numeric fields.", type="warning")
            return
 
        status.clear()
        run_btn.props("loading")
        run_btn.disable()
        try:
            r = await ng_run.io_bound(
                run_snr,
                aperture.value, temperature.value, diff_limit.value, mode.value,
                extime.value, template.value, is_model(),
                distance.value if is_model() else None,
            )
        except Exception as e:  # show errors in the page instead of crashing
            state["result"] = None
            dl_txt.disable()
            dl_html.disable()
            show_message(f"{type(e).__name__}: {e}", "error")
        else:
            state["result"] = r
            plot.figure = make_figure(r)
            plot.update()
            dl_txt.enable()
            dl_html.enable()
            if r["warning"]:
                show_message("WARNING: " + r["warning"], "warn")
            show_message(
                f"S/N Ratio computed for HWO with {r['pollux'].name} in mode "
                f"{r['pollux'].mode_name} for target {r['template']}.", "ok")
        finally:
            run_btn.props(remove="loading")
            run_btn.enable()
 
    def download_txt():
        r = state["result"]
        if r is not None:
            ui.download.content(snr_table_bytes(r), f"hwo_{r['template']}_{r['mode']}_snr.txt",
                                "text/plain")
 
    def download_html():
        r = state["result"]
        if r is not None:
            html = make_figure(r).to_html(include_plotlyjs="cdn").encode()
            ui.download.content(html, f"hwo_{r['template']}_{r['mode']}_snr.html", "text/html")
 
    run_btn.on_click(compute)
    dl_txt.on_click(download_txt)
    dl_html.on_click(download_html)
 
 
def main(host="127.0.0.1", port=8000, show=True, reload=False):
    ui.run(host=host, port=port, title="POLLUX ETC", favicon="🔭", show=show, reload=reload)
 
 
if __name__ in {"__main__", "__mp_main__"}:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--reload", action="store_true")
    args, _ = parser.parse_known_args()
    main(args.host, args.port, not args.no_browser, args.reload)