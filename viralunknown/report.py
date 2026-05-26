import yaml
import json
from pathlib import Path
from jinja2 import Template
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots





def make_kraken_pie(sample, taxa):
    if not taxa:
        return ""

    # Garder uniquement les espèces (rank S)
    species_only = {
        name: data for name, data in taxa.items()
        if data["rank"] == "S"
    }

    if not species_only:
        return ""

    sorted_taxa = sorted(species_only.items(), key=lambda x: x[1]["pct"], reverse=True)[:15]
    total_pct   = sum(t[1]["pct"] for t in sorted_taxa)
    other_pct   = max(0, 100 - total_pct)

    labels = [t[0] for t in sorted_taxa]
    values = [t[1]["pct"] for t in sorted_taxa]

    if other_pct > 0.1:
        labels.append("Others")
        values.append(other_pct)

    colors = px.colors.qualitative.Plotly + px.colors.qualitative.Set2

    fig = go.Figure(go.Pie(
        labels        = labels,
        values        = values,
        hole          = 0.4,
        textinfo      = "label+percent",
        hovertemplate = "<b>%{label}</b><br>%{value:.2f}%<extra></extra>",
        marker        = dict(colors=colors[:len(labels)]),
    ))
    fig.update_layout(
        title         = f"Kraken2 species — {sample}",
        height        = 500,
        showlegend    = True,
        legend        = dict(orientation="v", x=1.05),
        paper_bgcolor = "rgba(0,0,0,0)",
        font          = dict(size=12, color="#e2e8f0"),
    )
    return fig.to_html(full_html=False, include_plotlyjs=False)


def collect_stats(config):
    """Collect statistics from each pipeline step logs."""
    with open(config) as f:
        cfg = yaml.safe_load(f)

    samples = list(cfg["samples"].keys())
    stats   = {s: {} for s in samples}

    for sample in samples:

        # ── Raw reads (seqkit stats) ──────────────────────────────
        stat_file = Path(f"results/stats/{sample}_raw_stats.txt")
        if stat_file.exists():
            lines = stat_file.read_text().splitlines()
            if len(lines) >= 2:
                parts = lines[1].split()
                try:
                    stats[sample]["raw_reads"] = int(parts[3].replace(",", ""))
                    stats[sample]["raw_bases"]  = int(parts[4].replace(",", ""))
                    stats[sample]["avg_len"]    = float(parts[6].replace(",", ""))
                    stats[sample]["avg_qual"]   = float(parts[14].replace(",", ""))
                except (IndexError, ValueError):
                    pass

        # ── Dustmasker ────────────────────────────────────────────
        log = Path(f"logs/dustmasker/{sample}.log")
        if log.exists():
            for line in log.read_text().splitlines():
                if "complexité" in line or "complexity" in line.lower():
                    try:
                        stats[sample]["after_dust"] = int(line.split(":")[1].strip())
                    except (IndexError, ValueError):
                        pass

        # ── Kraken2 taxonomy ──────────────────────────────────────
        kreport = Path(f"results/kraken2/{sample}.kreport2")
        if kreport.exists():
            taxa = {}
            for line in kreport.read_text().splitlines():
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                try:
                    pct   = float(parts[0])
                    reads = int(parts[1])
                    rank  = parts[3]
                    name  = parts[5].strip()
                    if rank in ("F", "G", "S") and pct >= 0.1:
                        taxa[name] = {"pct": pct, "reads": reads, "rank": rank}
                except (ValueError, IndexError):
                    pass
            stats[sample]["taxa"] = taxa

            for line in kreport.read_text().splitlines():
                parts = line.strip().split("\t")
                if len(parts) >= 6 and parts[5].strip() == "root":
                    try:
                        stats[sample]["kraken_classified"] = int(parts[1])
                    except ValueError:
                        pass
                    break

        # ── Non-human reads ───────────────────────────────────────
        log = Path(f"logs/extract/{sample}.log")
        if log.exists():
            for line in log.read_text().splitlines():
                if "non humains" in line or "non-human" in line.lower():
                    try:
                        stats[sample]["non_human"] = int(line.split(":")[1].strip())
                    except (IndexError, ValueError):
                        pass

        # ── Unclassified reads ────────────────────────────────────
        log = Path(f"logs/unclassified/{sample}.log")
        if log.exists():
            for line in log.read_text().splitlines():
                if "non classifiés" in line or "unclassified" in line.lower():
                    try:
                        stats[sample]["unclassified"] = int(line.split(":")[1].strip())
                    except (IndexError, ValueError):
                        pass

        # ── Not mapped after Flye ─────────────────────────────────
        log = Path(f"logs/metaflye/{sample}_mapping.log")
        if log.exists():
            for line in log.read_text().splitlines():
                if "non mappés" in line or "unmapped" in line.lower():
                    try:
                        stats[sample]["not_mapped"] = int(line.split(":")[1].strip())
                    except (IndexError, ValueError):
                        pass

        # ── Quality filter (after Flye, before MMseqs2) ───────────
        log = Path(f"logs/quality_filter/{sample}.log")
        if log.exists():
            for line in log.read_text().splitlines():
                if "High quality" in line:
                    try:
                        stats[sample]["high_qual"] = int(line.split(":")[1].strip())
                    except (IndexError, ValueError):
                        pass
                if "Low quality" in line:
                    try:
                        stats[sample]["low_qual"] = int(line.split(":")[1].strip())
                    except (IndexError, ValueError):
                        pass

        # ── MMseqs2 clusters ──────────────────────────────────────
        summary = Path(f"results/mmseqs/{sample}/clusters_summary.tsv")
        if summary.exists():
            lines = summary.read_text().splitlines()
            stats[sample]["clusters"] = max(0, len(lines) - 1)

    return samples, stats, cfg


def make_funnel_chart(sample, stats):
    """Funnel chart showing the progression of reads."""
    steps = [
        ("Raw reads",      stats.get("raw_reads",   0)),
        ("Filter dustmasker",      stats.get("after_dust",   0)),
        ("Non humans",           stats.get("non_human",    0)),
        ("Unclassified",        stats.get("unclassified", 0)),
        ("Not mapped",         stats.get("not_mapped",   0)),
        ("Quality filter (Q≥10)",   stats.get("high_qual",    0)),
        ("Clusters MMseqs2",      stats.get("clusters",     0)),
        
    ]
    steps = [(k, v) for k, v in steps if v > 0]
    if not steps:
        return ""

    labels = [s[0] for s in steps]
    values = [s[1] for s in steps]

    fig = go.Figure(go.Funnel(
        y = labels,
        x = values,
        textinfo = "value+percent initial",
        marker = dict(color=[
            "#3bf67d","#3b82f6", "#6366f1", "#8b5cf6",
            "#a855f7", "#ec4899", "#f43f5e",
        ])
    ))
    fig.update_layout(
        title = f"Reads progression — {sample}",
        height = 400,
        margin = dict(l=20, r=20, t=40, b=20),
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor  = "rgba(0,0,0,0)",
        font = dict(size=13),
    )
    return fig.to_html(full_html=False, include_plotlyjs=False)


def make_bar_comparison(samples, stats):
    """Comparison chart of samples."""
    keys = ["high_qual", "after_dust", "non_human", "unclassified", "not_mapped"]
    labels = ["Q≥10", "Post-dust", "Non humans", "Unclassified", "Not mapped"]
    colors = ["#3b82f6", "#6366f1", "#8b5cf6", "#a855f7", "#f43f5e"]

    fig = go.Figure()
    for key, label, color in zip(keys, labels, colors):
        fig.add_trace(go.Bar(
            name = label,
            x    = samples,
            y    = [stats[s].get(key, 0) for s in samples],
            marker_color = color,
        ))

    fig.update_layout(
        barmode = "group",
        title   = "Step-by-step comparison of samples",
        height  = 450,
        xaxis_title = "Samples",
        yaxis_title = "Number of reads",
        legend  = dict(orientation="h", y=-0.2),
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor  = "rgba(0,0,0,0)",
        font = dict(size=12),
    )
    return fig.to_html(full_html=False, include_plotlyjs=False)


HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ViralUnknown — Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: 'Segoe UI', system-ui, sans-serif;
            background: #0f172a;
            color: #e2e8f0;
            line-height: 1.6;
        }
        header {
            background: linear-gradient(135deg, #1e40af, #7c3aed);
            padding: 2rem 3rem;
            border-bottom: 1px solid #334155;
        }
        header h1 { font-size: 2rem; font-weight: 700; }
        header p  { color: #bfdbfe; margin-top: 0.3rem; }
        .badge {
            display: inline-block;
            background: rgba(255,255,255,0.15);
            padding: 0.2rem 0.8rem;
            border-radius: 999px;
            font-size: 0.8rem;
            margin-top: 0.5rem;
        }
        main { max-width: 1400px; margin: 0 auto; padding: 2rem 3rem; }
        h2 {
            font-size: 1.3rem;
            font-weight: 600;
            color: #93c5fd;
            border-left: 4px solid #3b82f6;
            padding-left: 0.8rem;
            margin: 2rem 0 1rem;
        }
        h3 {
            font-size: 1.1rem;
            color: #c4b5fd;
            margin: 1.5rem 0 0.5rem;
        }
        .grid-4 {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
            margin-bottom: 1.5rem;
        }
        .card {
            background: #1e293b;
            border: 1px solid #334155;
            border-radius: 12px;
            padding: 1.2rem 1.5rem;
        }
        .card .label {
            font-size: 0.8rem;
            color: #94a3b8;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }
        .card .value {
            font-size: 1.8rem;
            font-weight: 700;
            color: #60a5fa;
            margin-top: 0.2rem;
        }
        .card .sub {
            font-size: 0.8rem;
            color: #64748b;
            margin-top: 0.2rem;
        }
        .chart-box {
            background: #1e293b;
            border: 1px solid #334155;
            border-radius: 12px;
            padding: 1rem;
            margin-bottom: 1.5rem;
        }
        .sample-section {
            background: #162032;
            border: 1px solid #334155;
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 2rem;
        }
        .sample-title {
            font-size: 1.2rem;
            font-weight: 700;
            color: #f0abfc;
            padding-bottom: 0.8rem;
            border-bottom: 1px solid #334155;
            margin-bottom: 1rem;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9rem;
        }
        th {
            background: #0f172a;
            padding: 0.6rem 1rem;
            text-align: left;
            color: #94a3b8;
            font-weight: 500;
        }
        td {
            padding: 0.6rem 1rem;
            border-top: 1px solid #334155;
        }
        tr:hover td { background: #243147; }
        .pct { color: #a78bfa; font-size: 0.85rem; }
        .two-col {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 1rem;
        }
        @media (max-width: 900px) {
            .two-col { grid-template-columns: 1fr; }
        }
        footer {
            text-align: center;
            padding: 2rem;
            color: #475569;
            font-size: 0.85rem;
            border-top: 1px solid #1e293b;
        }
    </style>
</head>
<body>
<header>
    <h1>🦠 ViralUnknown</h1>
    <p>Unknown virus detection report — Nanopore data</p>
    <span class="badge">{{ n_samples }} sample(s) • Generated on {{ date }}</span>
</header>

<main>

    <h2>Summary table</h2>
    <div class="chart-box">
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>Raw reads</th>
                <th>After dustmasker</th>
                <th>Non-human</th>
                <th>Unclassified</th>
                <th>Not assembled</th>
                <th>High quality</th>
                <th>Clusters</th>
            </tr>
        </thead>
        <tbody>
        {% for sample in samples %}
            <tr>
                <td>{{ sample }}</td>
                <td>{{ stats[sample].get('raw_reads', '—') }}</td>
                <td>
                    {{ stats[sample].get('after_dust', '—') }}
                    {% if stats[sample].get('raw_reads') and stats[sample].get('after_dust') %}
                    <span class="pct">({{ "%.1f"|format(stats[sample]['after_dust'] / stats[sample]['raw_reads'] * 100) }}%)</span>
                    {% endif %}
                </td>
                <td>{{ stats[sample].get('non_human', '—') }}</td>
                <td>{{ stats[sample].get('unclassified', '—') }}</td>
                <td>{{ stats[sample].get('not_mapped', '—') }}</td>
                <td>{{ stats[sample].get('high_qual', '—') }}</td>
                <td>{{ stats[sample].get('clusters', '—') }}</td>
            </tr>
        {% endfor %}
        </tbody>
    </table>
    </div>

    <h2>Per-sample details</h2>

    {% for sample in samples %}
    <div class="sample-section">
        <div class="sample-title">📊 {{ sample }}</div>

        <!-- Stats cards -->
        <div class="grid-4">
            <div class="card">
                <div class="label">Raw reads</div>
                <div class="value">{{ stats[sample].get('raw_reads', '—') }}</div>
            </div>
            <div class="card">
                <div class="label">Average length</div>
                <div class="value">{{ stats[sample].get('avg_len', '—') }}</div>
                <div class="sub">bases</div>
            </div>
            <div class="card">
                <div class="label">Q30</div>
                <div class="value">{{ stats[sample].get('avg_qual', '—') }}</div>
            </div>
            <div class="card">
                <div class="label">Final clusters</div>
                <div class="value">{{ stats[sample].get('clusters', '—') }}</div>
                <div class="sub">candidates to BLAST</div>
            </div>
        </div>

        <!-- Funnel + Pie side by side -->
        <div class="two-col">
            <div class="chart-box">{{ funnels[sample] }}</div>
            {% if kraken_pies[sample] %}
            <div class="chart-box">{{ kraken_pies[sample] }}</div>
            {% else %}
            <div class="chart-box" style="display:flex;align-items:center;justify-content:center;color:#475569;">
                No Kraken2 species data available
            </div>
            {% endif %}
        </div>

    </div>
    {% endfor %}

</main>

<footer>
    ViralUnknown v{{ version }} — Nanopore viral metagenomics pipeline
</footer>
</body>
</html>
"""

def generate_report(config, output):
    from datetime import datetime
    from . import __version__

    samples, stats, cfg = collect_stats(config)

    # Camemberts Kraken2
    kraken_pies = {
        s: make_kraken_pie(s, stats[s].get("taxa", {}))
        for s in samples
    }

    total_reads        = sum(s.get("raw_reads",    0) for s in stats.values())
    total_unclassified = sum(s.get("unclassified", 0) for s in stats.values())
    total_clusters     = sum(s.get("clusters",     0) for s in stats.values())
    pct_unclassified   = (
        f"{total_unclassified / total_reads * 100:.1f}"
        if total_reads > 0 else "0"
    )

    funnels   = {s: make_funnel_chart(s, stats[s]) for s in samples}
    bar_chart = make_bar_comparison(samples, stats)

    template = Template(HTML_TEMPLATE)
    html = template.render(
        samples            = samples,
        stats              = stats,
        n_samples          = len(samples),
        total_reads        = f"{total_reads:,}",
        total_unclassified = f"{total_unclassified:,}",
        total_clusters     = f"{total_clusters:,}",
        pct_unclassified   = pct_unclassified,
        bar_chart          = bar_chart,
        funnels            = funnels,
        kraken_pies        = kraken_pies,      # ← vérifie que cette ligne est là
        date               = datetime.now().strftime("%d/%m/%Y %H:%M"),
        version            = __version__,
    )

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(html)