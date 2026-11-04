"""SWAM profile editor (Textual), styled after an old-school MMO client.

Layout: orb status bar on top, main viewport (one tab at a time) with a stone
side-tab panel on the right, a parchment "chatbox" log below, footer keys last.
Only the user's profile YAML is edited; comments are preserved (ruamel round-trip).
"""

from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import Any

from .config import create_backup_for, default_profile_resource, prune_empty_profile_sections

try:
    from ruamel.yaml import YAML
    from textual.app import App, ComposeResult
    from textual.containers import Grid, Horizontal, Vertical, VerticalScroll
    from textual.theme import Theme
    from textual.widgets import Button, ContentSwitcher, Footer, Input, Label, RichLog, Select, Static, Switch
except ModuleNotFoundError as exc:  # pragma: no cover - depends on user environment.
    MISSING_TUI_DEPENDENCY = exc.name or "a required dependency"
else:
    MISSING_TUI_DEPENDENCY = None

# --------------------------------------------------------------------------
#  Field definitions
# --------------------------------------------------------------------------

# (config key, tile name, glyph, description)
STAGE_TILES = [
    ("run_short_reads", "Short reads", "≈", "fastp · host filter · KMA/DIAMOND · Nonpareil"),
    ("run_contigs", "Contigs", "§", "MEGAHIT · geNomad · MobMess · AMRFinder · MMseqs2"),
    ("run_mags", "MAGs", "♦", "MetaBAT2 · CheckM2 · GTDB-tk · METABOLIC · CoverM"),
]

MAG_SKIP_FIELDS = [
    ("skip_gtdbtk", "Skip GTDB-tk", "needs gtdbtk_db"),
    ("skip_metabolic", "Skip METABOLIC", ""),
    ("skip_checkm2", "Skip CheckM2", ""),
]

STAGE_INT_FIELDS = [
    ("genomad_splits", "geNomad DB splits", "raise to lower geNomad memory"),
    ("min_read_length", "Min read length (bp)", "fastp; ≥24 for Nonpareil"),
]

AMR_MODES = [("operational", "operational"), ("discovery", "discovery")]

PATH_FIELDS = [
    ("in_dir", "Input FASTQ dir", "paired-end *R1*/*_1 FASTQs"),
    ("out_dir", "Output dir", "analysis outputs + logs/"),
    ("dbs_dir", "Database dir", "empty = ~/.local/share/swam/dbs"),
    ("gtdbtk_db", "GTDB-tk data dir", "manual download; MAG taxonomy"),
    ("clean_reads_dir", "Clean reads dir", "only if Short reads is off"),
    ("contigs_dir", "Contigs dir", "only if Contigs is off ({sample}.contigs.fa)"),
]

# Keep in sync with common.smk sample discovery (and build_run_report.py SAMPLE_REGEX).
R1_GLOBS = ("*R1*.fastq*", "*_1.fastq*")

SNAKEMAKE_INT_FIELDS = [
    ("jobs", "Max submitted jobs"),
    ("cores", "Total cores"),
    ("restart-times", "Restart attempts"),
    ("latency-wait", "Latency wait (s)"),
]

SNAKEMAKE_BOOL_FIELDS = [
    ("use-conda", "Use conda"),
    ("keep-going", "Keep going after failures"),
    ("keep-incomplete", "Keep incomplete files"),
    ("rerun-incomplete", "Rerun incomplete jobs"),
    ("show-failed-logs", "Show failed logs"),
]

CONDA_FRONTENDS = [("conda", "conda"), ("mamba", "mamba")]

SLURM_DEFAULT_FIELDS = [
    ("slurm_account", "Slurm account", "required by most clusters"),
    ("slurm_partition", "Slurm partition", ""),
    ("mem_mb", "Default mem (MB)", "rules without an override"),
    ("runtime", "Default walltime (min)", ""),
]

RULE_ORDER = [
    "short_reads",
    "short_reads_argnet",
    "contigs",
    "genomad",
    "mobmess",
    "prodigal",
    "contig_amr",
    "contig_argnet",
    "mge_annotation",
    "mmseqs_taxonomy",
    "contig_abundance",
    "bin",
    "mag_prodigal",
    "mag_amr",
    "mag_argnet",
    "mag_mge",
    "mag_qc",
    "mag_taxonomy",
    "mag_metabolism",
    "mag_abundance",
    "init_mmseqs_db",
]

# Keep in sync with `localrules:` in workflow/Snakefile — these run on the head
# node, so Slurm resources don't apply; their YAML entries are left untouched.
LOCAL_RULES = {
    "initiate_dbs",
    "short_reads_summary",
    "argnet_short_reads_summary",
    "init_genomad",
    "init_amrfinder_db",
    "init_argnet",
    "init_checkm2_db",
    "init_mge_tool",
    "contig_summary",
    "argnet_contig_summary",
    "assembly_qa",
    "amr_risk_score",
    "mag_summary",
    "argnet_mag_summary",
}

# mem/runtime scale per retry attempt inside the rule files; only threads apply here.
RETRY_SCALED_RULES = {"contigs", "genomad", "bin"}

# --------------------------------------------------------------------------
#  Themes
# --------------------------------------------------------------------------

OSRS_THEME_VARS = {
    "stone": "#3e3529",
    "stone-dark": "#1e1a15",
    "stone-light": "#6f6455",
    "stone-shadow": "#0f0d0a",
    "gold": "#ffff00",
    "parchment": "#c8b28a",
    "parchment-ink": "#1a1206",
    "examine": "#00ffff",
}

CLASSIC_THEME_VARS = {
    "stone": "#111827",
    "stone-dark": "#020617",
    "stone-light": "#334155",
    "stone-shadow": "#000000",
    "gold": "#38bdf8",
    "parchment": "#0f172a",
    "parchment-ink": "#e5e7eb",
    "examine": "#94a3b8",
}

CHAT_COLORS = {
    "osrs": {"info": "#1a1206", "ok": "#006400", "warn": "#7a3d00", "error": "#b00000", "game": "#0000c0"},
    "classic": {"info": "#e5e7eb", "ok": "#22c55e", "warn": "#facc15", "error": "#ef4444", "game": "#38bdf8"},
}


def build_themes() -> dict[str, "Theme"]:
    return {
        "osrs": Theme(
            name="swam-osrs",
            primary="#ff981f",
            secondary="#ffff00",
            accent="#ff981f",
            foreground="#e8dcc0",
            background="#2b2620",
            surface="#3e3529",
            panel="#1e1a15",
            success="#00c800",
            warning="#ff981f",
            error="#ff3030",
            dark=True,
            variables=OSRS_THEME_VARS,
        ),
        "classic": Theme(
            name="swam-classic",
            primary="#38bdf8",
            secondary="#1d4ed8",
            accent="#38bdf8",
            foreground="#e5e7eb",
            background="#0f172a",
            surface="#111827",
            panel="#020617",
            success="#22c55e",
            warning="#facc15",
            error="#ef4444",
            dark=True,
            variables=CLASSIC_THEME_VARS,
        ),
    }


def require_tui_dependencies() -> None:
    if MISSING_TUI_DEPENDENCY is None:
        return
    missing = {"ruamel": "ruamel.yaml"}.get(MISSING_TUI_DEPENDENCY, MISSING_TUI_DEPENDENCY)
    raise RuntimeError(
        f"Missing dependency: {missing}\n\n"
        "Install the SWAM CLI dependencies, then run the command again:\n"
        "  conda install -c conda-forge textual ruamel.yaml\n"
        "or\n"
        "  python -m pip install -e .\n"
    )


# --------------------------------------------------------------------------
#  Helpers
# --------------------------------------------------------------------------


def load_yaml(path: Path) -> tuple[Any, Any]:
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.load(handle)
    return yaml, data


def load_yaml_text(text: str) -> Any:
    yaml = YAML()
    yaml.preserve_quotes = True
    return yaml.load(text)


def as_text(value: Any) -> str:
    return "" if value is None else str(value)


def count_samples(in_dir: str) -> int | None:
    """Number of R1 FASTQs in in_dir, or None when the directory doesn't exist."""
    if not in_dir:
        return None
    path = Path(in_dir).expanduser()
    if not path.is_dir():
        return None
    hits: set[str] = set()
    for pattern in R1_GLOBS:
        hits.update(glob.glob(os.path.join(str(path), pattern)))
    return len(hits)


def parse_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{label} must be a whole number (got '{value}').") from exc
    if parsed < 0:
        raise ValueError(f"{label} must be zero or greater.")
    return parsed


def parse_runtime(value: str, label: str) -> int | str:
    """Minutes as int when numeric; otherwise keep Slurm-style strings like '6h'."""
    return parse_int(value, label) if value.isdigit() else value


if MISSING_TUI_DEPENDENCY is None:

    class FieldRow(Horizontal):
        def __init__(self, label: str, control: Any, help_text: str = "") -> None:
            super().__init__(classes="field-row")
            self.label_text = label
            self.control = control
            self.help_text = help_text

        def compose(self) -> ComposeResult:
            yield Label(self.label_text, classes="field-label")
            yield self.control
            yield Static(self.help_text, classes="field-help")

    class Section(Vertical):
        def __init__(self, title: str, *children: Any, id: str | None = None) -> None:
            super().__init__(classes="section", id=id)
            self.border_title = title
            self.section_children = children

        def compose(self) -> ComposeResult:
            yield from self.section_children

    class StageTile(Vertical):
        """A 'prayer' tile: glows gold when its stage is switched on."""

        def __init__(self, key: str, name: str, glyph: str, desc: str, value: bool) -> None:
            super().__init__(classes="tile", id=f"tile__{key}")
            self.key = key
            self.tile_name = name
            self.glyph = glyph
            self.desc = desc
            self.value = value

        def compose(self) -> ComposeResult:
            yield Static(self.glyph, classes="tile-glyph")
            yield Static(self.tile_name, classes="tile-name")
            yield Static(self.desc, classes="tile-desc")
            yield Switch(value=self.value, id=f"config__{self.key}")

        def on_mount(self) -> None:
            self.set_class(self.value, "active")

    class ConfigEditorApp(App[None]):
        TITLE = "SWAM"
        CSS_PATH = "tui.tcss"

        TABS = [
            ("stages", "†", "Stages"),
            ("paths", "✎", "Paths"),
            ("slurm", "♜", "Slurm"),
            ("snakemake", "≡", "Snakemake"),
        ]

        BINDINGS = [
            ("ctrl+s", "save", "Save"),
            ("ctrl+r", "reload", "Reload"),
            ("1", "show_tab('stages')", "Stages"),
            ("2", "show_tab('paths')", "Paths"),
            ("3", "show_tab('slurm')", "Slurm"),
            ("4", "show_tab('snakemake')", "Snakemake"),
            ("q", "quit_app", "Quit"),
            ("ctrl+q", "quit_app", "Quit"),
        ]

        def __init__(self, config_path: Path, profile: str, theme_name: str = "osrs") -> None:
            super().__init__()
            self.config_path = config_path
            self.profile = profile
            self.theme_key = theme_name
            self.dirty = False
            self.quit_armed = False
            # Themes must exist before the stylesheet (which uses their variables) is parsed.
            themes = build_themes()
            for theme in themes.values():
                self.register_theme(theme)
            self.theme = themes[theme_name].name
            self._load()

        # ---------------------------------------------------------- data ---

        def _load(self) -> None:
            self.yaml, self.data = load_yaml(self.config_path)
            self.default_data = load_yaml_text(default_profile_resource(self.profile).read_text(encoding="utf-8"))
            self.data.setdefault("config", {})
            self.is_slurm = str(self.data.get("executor", "")) == "slurm"
            self.rule_names = self._collect_rule_names()

        def cfg(self, key: str, fallback: Any = "") -> Any:
            value = self.data["config"].get(key)
            if value is None:
                value = (self.default_data.get("config") or {}).get(key, fallback)
            return value

        def _collect_rule_names(self) -> list[str]:
            found = set(RULE_ORDER)
            for source in (self.data, self.default_data):
                found.update((source.get("set-threads") or {}).keys())
                found.update((source.get("set-resources") or {}).keys())
            found -= LOCAL_RULES
            return [r for r in RULE_ORDER if r in found] + sorted(found - set(RULE_ORDER))

        def _section_value(self, section: str, rule: str, field: str | None = None) -> str:
            for source in (self.data, self.default_data):
                block = (source.get(section) or {}).get(rule)
                if block is None:
                    continue
                if field is None:
                    return as_text(block)
                if field in block:
                    return as_text(block[field])
            return ""

        def rule_group(self, rule: str) -> str | None:
            return (self.data.get("groups") or {}).get(rule) or (self.default_data.get("groups") or {}).get(rule)

        # ------------------------------------------------------- compose ---

        def compose(self) -> ComposeResult:
            with Horizontal(id="orbs"):
                yield Static(id="orb-samples", classes="orb")
                yield Static(id="orb-stages", classes="orb")
                yield Static(f"● Profile\n{self.profile}", id="orb-profile", classes="orb")
                yield Static(id="orb-amr", classes="orb")
                yield Static("SWAM", id="brand")
            with Horizontal(id="client"):
                with ContentSwitcher(initial="tab-stages", id="viewport"):
                    with VerticalScroll(id="tab-stages", classes="tab-body"):
                        yield from self.stages_tab()
                    with VerticalScroll(id="tab-paths", classes="tab-body"):
                        yield from self.paths_tab()
                    with VerticalScroll(id="tab-slurm", classes="tab-body"):
                        yield from self.slurm_tab()
                    with VerticalScroll(id="tab-snakemake", classes="tab-body"):
                        yield from self.snakemake_tab()
                with Vertical(id="sidepanel"):
                    yield Static("Settings", id="side-title")
                    with Grid(id="side-tabs"):
                        for key, glyph, name in self.TABS:
                            yield Button(f"{glyph}\n{name}", id=f"sidetab-{key}", classes="side-tab")
                    yield Static(self.side_info(), id="side-info")
                    yield Button("Save", id="save", classes="stone-button")
                    yield Button("Reload", id="reload", classes="stone-button")
            yield RichLog(id="chatbox", markup=True, wrap=True, max_lines=500)
            yield Footer()

        def side_info(self) -> str:
            return (
                "[b]Profile file[/b]\n"
                f"{self.config_path}\n\n"
                "[b]Run with[/b]\n"
                f"swam run --profile {self.profile}"
            )

        def stages_tab(self) -> list[Any]:
            tiles = Grid(
                *(StageTile(key, name, glyph, desc, bool(self.cfg(key, True))) for key, name, glyph, desc in STAGE_TILES),
                id="tiles",
            )
            amr = FieldRow(
                "AMR mode",
                Select(AMR_MODES, value=str(self.cfg("amr_mode", "operational")), id="config__amr_mode", allow_blank=False),
                "discovery adds ARGNet summaries",
            )
            skips = [
                FieldRow(label, Switch(value=bool(self.cfg(key, False)), id=f"config__{key}"), hint)
                for key, label, hint in MAG_SKIP_FIELDS
            ]
            ints = [
                FieldRow(label, Input(value=as_text(self.cfg(key)), id=f"config__{key}", type="integer"), hint)
                for key, label, hint in STAGE_INT_FIELDS
            ]
            return [
                Section("Pipeline stages", tiles, Static("", id="stage-hint", classes="hint")),
                Section("Annotation", amr, *ints),
                Section("MAG options", *skips),
            ]

        def paths_tab(self) -> list[Any]:
            rows = [
                FieldRow(label, Input(value=as_text(self.cfg(key)), id=f"config__{key}", placeholder="(not set)"), hint)
                for key, label, hint in PATH_FIELDS
            ]
            return [Section("Inputs & outputs", *rows)]

        def slurm_tab(self) -> list[Any]:
            if not self.is_slurm:
                return [
                    Section(
                        "Slurm",
                        Static(
                            "Slurm is not used by the local profile.\n"
                            "Edit it with:  swam config --profile slurm",
                            classes="hint",
                        ),
                    )
                ]
            defaults = self.data.setdefault("default-resources", {})
            default_rows = [
                FieldRow(label, Input(value=as_text(defaults.get(key)), id=f"default__{key}"), hint)
                for key, label, hint in SLURM_DEFAULT_FIELDS
            ]
            header = Horizontal(
                Static("Rule", classes="rule-name"),
                Static("Threads", classes="rule-col"),
                Static("Mem MB", classes="rule-col"),
                Static("Minutes", classes="rule-col"),
                Static("Group", classes="rule-col"),
                Static("Notes", classes="rule-note"),
                classes="rule-row rule-header",
            )
            rows: list[Any] = [
                Static(
                    "Group size batches N per-sample jobs into one Slurm job (run sequentially), "
                    "so walltime effectively multiplies. Blank = use defaults.",
                    classes="hint",
                ),
                header,
            ]
            for rule in self.rule_names:
                retry = rule in RETRY_SCALED_RULES
                group = self.rule_group(rule)
                group_value = self._section_value("group-components", group) if group else ""
                note = "mem/time retry-scaled in rules" if retry else (f"group '{group}'" if group else "")
                rows.append(
                    Horizontal(
                        Static(rule, classes="rule-name"),
                        Input(value=self._section_value("set-threads", rule), id=f"threads__{rule}", classes="rule-input"),
                        Input(
                            value="" if retry else self._section_value("set-resources", rule, "mem_mb"),
                            id=f"mem__{rule}",
                            classes="rule-input",
                            disabled=retry,
                        ),
                        Input(
                            value="" if retry else self._section_value("set-resources", rule, "runtime"),
                            id=f"runtime__{rule}",
                            classes="rule-input",
                            disabled=retry,
                        ),
                        Input(value=group_value, id=f"group__{rule}", classes="rule-input", disabled=group is None),
                        Static(note, classes="rule-note"),
                        classes="rule-row",
                    )
                )
            return [Section("Default resources", *default_rows), Section("Per-rule resources", *rows)]

        def snakemake_tab(self) -> list[Any]:
            rows: list[Any] = [
                FieldRow("Executor", Static(as_text(self.data.get("executor")), classes="readonly"), "set by profile"),
                FieldRow("Scheduler", Static("greedy", classes="readonly"), "always greedy (no cbc solver)"),
            ]
            for key, label in SNAKEMAKE_INT_FIELDS:
                if key in self.data:
                    rows.append(FieldRow(label, Input(value=as_text(self.data.get(key)), id=f"top__{key}", type="integer")))
            if "conda-frontend" in self.data:
                rows.append(
                    FieldRow(
                        "Conda frontend",
                        Select(CONDA_FRONTENDS, value=str(self.data.get("conda-frontend")), id="top__conda-frontend", allow_blank=False),
                    )
                )
            switches = [
                FieldRow(label, Switch(value=bool(self.data.get(key)), id=f"top__{key}"))
                for key, label in SNAKEMAKE_BOOL_FIELDS
                if key in self.data
            ]
            sections = [Section("Snakemake", *rows), Section("Behaviour", *switches)]
            if not self.is_slurm:
                pool = self.data.setdefault("resources", {})
                sections.append(
                    Section(
                        "Local resource pool",
                        FieldRow(
                            "Memory pool (MB)",
                            Input(value=as_text(pool.get("mem_mb")), id="resources__mem_mb", type="integer"),
                            "summed mem_mb of running jobs never exceeds this",
                        ),
                    )
                )
            return sections

        # ------------------------------------------------------ lifecycle ---

        def on_mount(self) -> None:
            self.query_one("#chatbox").border_title = "Chat"
            self.query_one("#viewport").border_title = "Stages"
            self.chat("Welcome to SWAM.", "game")
            self.chat(f"Editing the [b]{self.profile}[/b] profile.")
            self.mark_active_tab("stages")
            self.refresh_status()
            self.call_after_refresh(self._arm_dirty_tracking)

        def _arm_dirty_tracking(self) -> None:
            self.tracking = True

        tracking = False

        def chat(self, message: str, level: str = "info") -> None:
            color = CHAT_COLORS[self.theme_key][level]
            self.query_one("#chatbox", RichLog).write(f"[{color}]{message}[/]")

        # --------------------------------------------------------- events ---

        def on_button_pressed(self, event: Button.Pressed) -> None:
            button_id = event.button.id or ""
            if button_id.startswith("sidetab-"):
                self.action_show_tab(button_id.removeprefix("sidetab-"))
            elif button_id == "save":
                self.action_save()
            elif button_id == "reload":
                self.action_reload()

        def on_switch_changed(self, event: Switch.Changed) -> None:
            switch_id = event.switch.id or ""
            if switch_id.startswith("config__run_"):
                key = switch_id.removeprefix("config__")
                self.query_one(f"#tile__{key}").set_class(event.value, "active")
            self._changed()

        def on_input_changed(self, event: Input.Changed) -> None:
            self._changed()

        def on_select_changed(self, event: Select.Changed) -> None:
            self._changed()

        def _changed(self) -> None:
            if self.tracking:
                self.dirty = True
                self.quit_armed = False
            self.refresh_status()

        # -------------------------------------------------------- status ---

        def _input(self, widget_id: str) -> str:
            return self.query_one(f"#{widget_id}", Input).value.strip()

        def _switch(self, widget_id: str) -> bool:
            return bool(self.query_one(f"#{widget_id}", Switch).value)

        def _select(self, widget_id: str) -> str:
            return str(self.query_one(f"#{widget_id}", Select).value)

        def stage_problems(self) -> list[str]:
            """Mirror the input checks in common.smk."""
            sr, ctg, mags = (self._switch(f"config__{key}") for key, *_ in STAGE_TILES)
            problems = []
            if not (sr or ctg or mags):
                problems.append("No stages enabled — nothing will run.")
            if not sr and (ctg or mags) and not self._input("config__clean_reads_dir"):
                problems.append("Short reads off: set 'Clean reads dir' on the Paths tab.")
            if not ctg and mags and not self._input("config__contigs_dir"):
                problems.append("Contigs off but MAGs on: set 'Contigs dir' on the Paths tab.")
            if mags and not self._switch("config__skip_gtdbtk") and not self._input("config__gtdbtk_db"):
                problems.append("GTDB-tk enabled without 'GTDB-tk data dir' — set it or skip GTDB-tk.")
            return problems

        def refresh_status(self) -> None:
            in_dir = self._input("config__in_dir")
            samples = count_samples(in_dir)
            orb = self.query_one("#orb-samples", Static)
            if samples is None:
                orb.update(f"● Samples\n{'not set' if not in_dir else 'no dir'}")
                orb.set_class(True, "orb-low")
            else:
                orb.update(f"● Samples\n{samples}")
                orb.set_class(samples == 0, "orb-low")
            enabled = sum(self._switch(f"config__{key}") for key, *_ in STAGE_TILES)
            stages_orb = self.query_one("#orb-stages", Static)
            stages_orb.update(f"● Stages\n{enabled}/3")
            stages_orb.set_class(enabled == 0, "orb-low")
            self.query_one("#orb-amr", Static).update(f"● AMR\n{self._select('config__amr_mode')}")
            problems = self.stage_problems()
            hint = self.query_one("#stage-hint", Static)
            hint.update("\n".join(f"! {p}" for p in problems) if problems else "All stage inputs look consistent.")
            hint.set_class(bool(problems), "hint-warn")

        # -------------------------------------------------------- actions ---

        def mark_active_tab(self, tab: str) -> None:
            for key, _glyph, _name in self.TABS:
                self.query_one(f"#sidetab-{key}").set_class(key == tab, "active")

        def action_show_tab(self, tab: str) -> None:
            self.query_one("#viewport", ContentSwitcher).current = f"tab-{tab}"
            self.query_one("#viewport").border_title = dict((k, n) for k, _g, n in self.TABS)[tab]
            self.mark_active_tab(tab)

        def action_reload(self) -> None:
            self._load()
            self.dirty = False
            self.tracking = False
            self.refresh(recompose=True)
            self.call_after_refresh(self._after_reload)

        def _after_reload(self) -> None:
            self.query_one("#chatbox").border_title = "Chat"
            self.mark_active_tab("stages")
            self.refresh_status()
            self.chat("Reloaded profile from disk.", "game")
            self.tracking = True

        def action_quit_app(self) -> None:
            if self.dirty and not self.quit_armed:
                self.quit_armed = True
                self.chat("You have unsaved changes. Press Q again to quit without saving, or Ctrl+S to save.", "warn")
                return
            self.exit()

        def action_save(self) -> None:
            try:
                self.apply_form_to_data()
            except ValueError as exc:
                self.chat(f"Not saved: {exc}", "error")
                return
            prune_empty_profile_sections(self.data)
            backup_path = create_backup_for(self.config_path, self.profile)
            with self.config_path.open("w", encoding="utf-8") as handle:
                self.yaml.dump(self.data, handle)
            self.dirty = False
            self.quit_armed = False
            self.chat(f"Saved {self.config_path.name} (backup: {backup_path}).", "ok")
            for warning in self.save_warnings():
                self.chat(warning, "warn")
            self.chat(f"Run with: swam run --profile {self.profile}", "game")

        def save_warnings(self) -> list[str]:
            warnings = list(self.stage_problems())
            for key, label, _hint in PATH_FIELDS:
                value = self._input(f"config__{key}")
                if value and key != "out_dir" and not Path(value).expanduser().exists():
                    warnings.append(f"{label} does not exist here: {value}")
            if count_samples(self._input("config__in_dir")) == 0:
                warnings.append("No *R1*/*_1 FASTQs found in the input dir.")
            if self.is_slurm and not self._input("default__slurm_account"):
                warnings.append("Slurm account is empty; many clusters reject jobs without one.")
            return warnings

        def apply_form_to_data(self) -> None:
            """Validate everything first, then write into the ruamel tree."""
            config = self.data["config"]
            updates: list[tuple[Any, str, Any]] = []

            for key, *_ in STAGE_TILES:
                updates.append((config, key, self._switch(f"config__{key}")))
            for key, *_ in MAG_SKIP_FIELDS:
                updates.append((config, key, self._switch(f"config__{key}")))
            updates.append((config, "amr_mode", self._select("config__amr_mode")))
            for key, label, _hint in STAGE_INT_FIELDS:
                updates.append((config, key, parse_int(self._input(f"config__{key}"), label)))
            for key, *_ in PATH_FIELDS:
                updates.append((config, key, self._input(f"config__{key}")))

            for key, label in SNAKEMAKE_INT_FIELDS:
                if key in self.data:
                    updates.append((self.data, key, parse_int(self._input(f"top__{key}"), label)))
            if "conda-frontend" in self.data:
                updates.append((self.data, "conda-frontend", self._select("top__conda-frontend")))
            for key, _label in SNAKEMAKE_BOOL_FIELDS:
                if key in self.data:
                    updates.append((self.data, key, self._switch(f"top__{key}")))

            if not self.is_slurm:
                pool = self.data.setdefault("resources", {})
                updates.append((pool, "mem_mb", parse_int(self._input("resources__mem_mb"), "Memory pool")))
            else:
                updates.extend(self._slurm_updates())

            for target, key, value in updates:
                if value is None:
                    target.pop(key, None)
                else:
                    target[key] = value

        def _slurm_updates(self) -> list[tuple[Any, str, Any]]:
            updates: list[tuple[Any, str, Any]] = []
            defaults = self.data.setdefault("default-resources", {})
            updates.append((defaults, "slurm_account", self._input("default__slurm_account")))
            updates.append((defaults, "slurm_partition", self._input("default__slurm_partition")))
            updates.append((defaults, "mem_mb", parse_int(self._input("default__mem_mb"), "Default mem")))
            updates.append((defaults, "runtime", parse_runtime(self._input("default__runtime"), "Default walltime")))

            set_threads = self.data.setdefault("set-threads", {})
            set_resources = self.data.setdefault("set-resources", {})
            group_components = self.data.setdefault("group-components", {})
            for rule in self.rule_names:
                threads = self._input(f"threads__{rule}")
                updates.append((set_threads, rule, parse_int(threads, f"{rule} threads") if threads else None))
                if rule not in RETRY_SCALED_RULES:
                    mem = self._input(f"mem__{rule}")
                    runtime = self._input(f"runtime__{rule}")
                    if mem or runtime:
                        block = set_resources.get(rule)
                        if block is None:
                            block = {}
                            updates.append((set_resources, rule, block))
                        updates.append((block, "mem_mb", parse_int(mem, f"{rule} mem") if mem else None))
                        updates.append((block, "runtime", parse_runtime(runtime, f"{rule} runtime") if runtime else None))
                    else:
                        updates.append((set_resources, rule, None))
                group = self.rule_group(rule)
                if group:
                    size = self._input(f"group__{rule}")
                    updates.append((group_components, group, parse_int(size, f"{rule} group size") if size else None))
            return updates


def launch_config_tui(config_path: Path, profile: str, theme_name: str = "osrs") -> int:
    require_tui_dependencies()
    ConfigEditorApp(config_path=config_path, profile=profile, theme_name=theme_name).run()
    return 0
