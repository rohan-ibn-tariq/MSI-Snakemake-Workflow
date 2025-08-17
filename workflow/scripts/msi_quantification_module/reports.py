"""
MSI Quantification Module - Reporting Functions

HTML report generation and data visualization.
Creates interactive charts and publication-ready HTML reports.
"""

import json
from datetime import datetime
import altair as alt
import polars as pl

alt.data_transformers.disable_max_rows()


def generate_msi_html_report(
    regional_results, af_evolution_results, msi_data, output_path
):
    """
    Generate MSI quantification HTML report.
    """

    overview_chart = create_overview_chart(msi_data)
    simple_indel_chart = create_simple_indel_chart(msi_data)
    indel_histogram = create_indel_histogram(msi_data)
    af_distribution = create_af_distribution_chart(msi_data)
    motif_breakdown_chart = create_motif_breakdown_chart(msi_data)

    # MSI Regional Analysis Charts & Summary Table
    regional_uncertainty_chart = create_msi_score_uncertainty_chart(regional_results)
    regional_breakdown_chart = create_region_count_breakdown_chart(regional_results)

    scope = msi_data["data_scope"]
    overview = msi_data["overview"]
    indels = msi_data["indels"]
    af = msi_data["allele_frequencies"]

    timestamp = datetime.now().strftime("%B %d, %Y")

    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MSI Quantification Report</title>
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
    <style>
        body {{
            font-family: "Times New Roman", Times, serif;
            background-color: #fefdfb;
            color: #1a1a1a;
            line-height: 1.6;
            max-width: 900px;
            margin: 0 auto;
            padding: 40px 20px;
        }}
        h1 {{
            text-align: center;
            font-size: 24px;
            font-weight: normal;
            border-bottom: 2px solid #1a1a1a;
            padding-bottom: 10px;
            margin-bottom: 30px;
        }}
        h2 {{
            font-size: 18px;
            font-weight: normal;
            margin-top: 40px;
            margin-bottom: 15px;
            # border-bottom: 1px solid #666;
        }}
        .header-info {{
            text-align: center;
            font-style: italic;
            margin-bottom: 40px;
            color: #666;
        }}
        .metrics {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .metric {{
            text-align: center;
            padding: 20px;
            border: 1px solid #ccc;
            background-color: #fbfaf8;
        }}
        .metric-value {{
            font-size: 32px;
            font-weight: bold;
            color: #8b4513;
            display: block;
        }}
        .metric-label {{
            font-size: 14px;
            color: #666;
            margin-top: 5px;
        }}
        .charts {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 30px;
            margin: 40px 0;
        }}
        .chart-container {{
            text-align: center;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-family: "Courier New", monospace;
            font-size: 12px;
        }}
        th, td {{
            border: 1px solid #999;
            padding: 8px;
            text-align: left;
        }}
        th {{
            background-color: #f5f4f2;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background-color: #fbfaf8;
        }}
        .note {{
            margin: 20px 0;
            padding: 15px;
            background-color: #f9f9f9;
            border-left: 4px solid #8b4513;
            font-style: italic;
            color: #666;
        }}
        .regional-dashboard {{
            margin: 40px 0;
            padding: 20px;
            background-color: #f9f9f9;
            border: 1px solid #ddd;
        }}
        
        .regional-charts {{
            display: grid;
            grid-template-columns: 1fr 1fr 1fr;
            gap: 20px;
            margin: 30px 0;
        }}
        
        .regional-chart-container {{
            text-align: center;
            background-color: white;
            padding: 15px;
            border: 1px solid #ccc;
        }}

    </style>
</head>
<body>
    <h1>MSI Quantification Report</h1>

    <div class="header-info">
        Microsatellite Instability Quantification Analysis<br>
        Generated {timestamp}
    </div>

    <!-- Regional MSI Analysis Dashboard -->
    <div class="regional-dashboard">
        <h2 style="color: #8b4513; font-weight: bold; font-size: 20px; margin: 40px 0 25px 0; text-align: left;">
            Complete Genome-Wide Analysis
        </h2>
        
        <div class="regional-charts">
            <div class="regional-chart-container">
                <div id="regional-uncertainty-chart"></div>
            </div>
            <div class="regional-chart-container">
                <div id="regional-breakdown-chart"></div>
            </div>
            <div class="regional-chart-container">
                <div id="regional-status-chart"></div>
            </div>
        </div>
        
        <!-- Regional analysis Summary Table -->
        <div class="regional-table-placeholder">
            <p><em>Regional analysis table coming next...</em></p>
        </div>
    </div>

    <div class="metrics">
        <div class="metric">
            <span class="metric-value">{overview['msi_rate']:.1f}%</span>
            <div class="metric-label">MSI Rate</div>
        </div>
        <div class="metric">
            <span class="metric-value">{overview['msi_regions']:,}</span>
            <div class="metric-label">MSI Unstable Regions</div>
        </div>
        <div class="metric">
            <span class="metric-value">{overview['non_msi_regions_with_indels']:,}</span>
            <div class="metric-label">Uncertain Regions</div>
        </div>
        <div class="metric">
            <span class="metric-value">{overview['stable_regions']:,}</span>
            <div class="metric-label">MSI Stable Regions</div>
        </div>
        <div class="metric">
            <span class="metric-value">{scope['total_ms_regions']:,}</span>
            <div class="metric-label">Total MS Regions</div>
        </div>
        <div class="metric">
            <span class="metric-value">{indels['insertion_count'] + indels['deletion_count']:,}</span>
            <div class="metric-label">Perfect MSI Indels</div>
        </div>
        <div class="metric">
            <span class="metric-value">{af['mean']:.2f}</span>
            <div class="metric-label">Mean AF_max (non-zero, non-N/A)</div>
        </div>
    </div>

    <div class="note">
        <strong>Analysis Scope:</strong> {scope['analysis_note']} 
        MSI stability refers to microsatellite-specific analysis only. Other variant types (SNVs, structural variants) were not analyzed.
    </div>

    <div class="charts">
        <div class="chart-container">
            <div id="overview-chart"></div>
        </div>
        <div class="chart-container">
            <div id="simple-indel-chart"></div>
        </div>
        <div class="chart-container">
            <div id="indel-histogram"></div>
        </div>
        <div class="chart-container">
            <div id="af-chart"></div>
        </div>
        <div class="chart-container">
            <div id="motif-chart"></div>
        </div>
    </div>




    <h2>Allele Frequency Analysis</h2>
    <table>
        <thead>
            <tr>
                <th>AF_max Category</th>
                <th>Count</th>
                <th>Percentage</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td>Variants with AF_max > 0</td>
                <td>{af['count']:,}</td>
                <td>{(af['count']/scope['perfect_msi_variants_analyzed']*100 if scope['perfect_msi_variants_analyzed'] > 0 else 0.0):.1f}%</td>
            </tr>
            <tr>
                <td>Variants with AF_max = 0</td>
                <td>{af['count_with_zero_af']:,}</td>
                <td>{(af['count_with_zero_af']/scope['perfect_msi_variants_analyzed']*100 if scope['perfect_msi_variants_analyzed'] > 0 else 0.0):.1f}%</td>
            </tr>
            <tr>
                <td>Variants with AF_max = N/A</td>
                <td>{af['count_with_na_af']:,}</td>
                <td>{(af['count_with_na_af']/scope['perfect_msi_variants_analyzed']*100 if scope['perfect_msi_variants_analyzed'] > 0 else 0.0):.1f}%</td>
            </tr>
            <tr>
                <td><strong>Total Perfect MSI Variants</strong></td>
                <td><strong>{scope['perfect_msi_variants_analyzed']:,}</strong></td>
                <td><strong>100.0%</strong></td>
            </tr>
        </tbody>
    </table>
    <div class="note">
        {af['mean_note']}
    </div>

    <h2>Motif Type Breakdown</h2>
    <table>
        <thead>
            <tr>
                <th>Motif Type</th>
                <th>Total Regions</th>
                <th>MSI Unstable</th>
                <th>MSI Rate</th>
                <th>Uncertain (N/A indels)</th>
                <th>MSI Stable</th>
            </tr>
        </thead>
        <tbody>"""

    for motif_type, data in sorted(msi_data["motif_breakdown"].items()):
        msi_rate = (data["msi"] / data["total"] * 100) if data["total"] > 0 else 0
        html_content += f"""
            <tr>
                <td>{motif_type}</td>
                <td>{data['total']:,}</td>
                <td>{data['msi']:,}</td>
                <td>{msi_rate:.1f}%</td>
                <td>{data['non_msi_with_indels']:,}</td>
                <td>{data['stable']:,}</td>
            </tr>"""

    html_content += f"""
        </tbody>
    </table>
    <div class="note">
        MSI Stable = no MSI indels detected in region. Uncertain = regions with N/A indels that may represent complex MSI events.
    </div>

    <script>
        vegaEmbed('#overview-chart', {overview_chart.to_json()}, {{actions: false}});
        vegaEmbed('#simple-indel-chart', {simple_indel_chart.to_json()}, {{actions: false}});
        vegaEmbed('#indel-histogram', {indel_histogram.to_json()}, {{actions: false}});
        vegaEmbed('#af-chart', {af_distribution.to_json()}, {{actions: false}});
        vegaEmbed('#motif-chart', {motif_breakdown_chart.to_json()}, {{actions: false}});
        
        // Regional MSI Analysis Charts
        vegaEmbed('#regional-uncertainty-chart', {regional_uncertainty_chart.to_json()}, {{actions: false}});
        vegaEmbed('#regional-breakdown-chart', {regional_breakdown_chart.to_json()}, {{actions: false}});  // ADD THIS LINE

    </script>
</body>
</html>
    """

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"[MSI-ANALYSIS] HTML report saved to: {output_path}")


def create_overview_chart(msi_data):
    """Create overview pie chart of region classification"""
    overview = msi_data["overview"]

    data = [
        {
            "Category": "MSI Unstable",
            "Count": overview["msi_regions"],
            "Color": "#8b4513",
        },
        {
            "Category": "Uncertain",
            "Count": overview["non_msi_regions_with_indels"],
            "Color": "#daa520",
        },
        {
            "Category": "MSI Stable",
            "Count": overview["stable_regions"],
            "Color": "#2f4f4f",
        },
    ]

    chart = (
        alt.Chart(pl.DataFrame(data))
        .mark_arc()
        .encode(
            theta=alt.Theta(field="Count", type="quantitative"),
            color=alt.Color(field="Color", type="nominal", scale=None),
            tooltip=[
                alt.Tooltip("Category:N", title="Category"),
                alt.Tooltip("Count:Q", title="Count", format=","),
            ],
        )
        .properties(width=200, height=200, title="MSI Region Classification")
    )

    return chart


def create_simple_indel_chart(msi_data):
    """Create simple histogram of insertion vs deletion counts with percentages"""
    indels = msi_data["indels"]
    total = indels["insertion_count"] + indels["deletion_count"]

    data = [
        {
            "Type": "Insertions",
            "Count": indels["insertion_count"],
            "Percentage": (indels["insertion_count"] / total * 100) if total > 0 else 0,
            "Label": (
                f"{indels['insertion_count']:,} ({(indels['insertion_count']/total*100):.1f}%)"
                if total > 0
                else "0"
            ),
            "Color": "#8b4513",
        },
        {
            "Type": "Deletions",
            "Count": indels["deletion_count"],
            "Percentage": (indels["deletion_count"] / total * 100) if total > 0 else 0,
            "Label": (
                f"{indels['deletion_count']:,} ({(indels['deletion_count']/total*100):.1f}%)"
                if total > 0
                else "0"
            ),
            "Color": "#2f4f4f",
        },
    ]

    bars = (
        alt.Chart(pl.DataFrame(data))
        .mark_bar(width=60)
        .encode(
            x=alt.X("Type:N", title="Indel Type"),
            y=alt.Y("Count:Q", title="Count"),
            color=alt.Color(field="Color", type="nominal", scale=None),
            tooltip=[
                alt.Tooltip("Type:N", title="Type"),
                alt.Tooltip("Count:Q", title="Count", format=","),
                alt.Tooltip("Percentage:Q", title="Percentage", format=".1f"),
            ],
        )
    )

    text = (
        alt.Chart(pl.DataFrame(data))
        .mark_text(align="center", baseline="bottom", dy=-5, fontSize=9)
        .encode(x=alt.X("Type:N"), y=alt.Y("Count:Q"), text=alt.Text("Label:N"))
    )

    chart = (bars + text).properties(
        width=200, height=150, title="Perfect MSI Indel Counts"
    )

    return chart


def create_indel_histogram(msi_data):
    """Create histogram of indel sizes with better hover and bins"""
    indels = msi_data["indels"]

    data = []
    for size in indels["insertion_sizes"]:
        data.append({"Size": size, "Type": "Insertion"})
    for size in indels["deletion_sizes"]:
        data.append({"Size": size, "Type": "Deletion"})

    chart = (
        alt.Chart(pl.DataFrame(data))
        .mark_bar(opacity=0.7, stroke="white", strokeWidth=0.5)
        .encode(
            x=alt.X("Size:Q", bin=alt.Bin(maxbins=15), title="Indel Size (bp)"),
            y=alt.Y("count():Q", title="Count"),
            color=alt.Color("Type:N", scale=alt.Scale(range=["#8b4513", "#2f4f4f"])),
            tooltip=[
                alt.Tooltip("Size:Q", title="Size Range (bp)", bin=alt.Bin(maxbins=15)),
                alt.Tooltip("count():Q", title="Count"),
                alt.Tooltip("Type:N", title="Type"),
            ],
        )
        .properties(width=200, height=150, title="MSI Indel Size Distribution")
    )

    return chart


def create_af_distribution_chart(msi_data):
    """Create AF_max distribution histogram (excluding 0 and N/A)"""
    af_values = msi_data["allele_frequencies"]["values"]

    if not af_values:
        return (
            alt.Chart(pl.DataFrame([{"x": 0, "y": 0}]))
            .mark_text(text="No AF data", size=14)
            .encode(x="x:Q", y="y:Q")
            .properties(width=200, height=150, title="AF_max Distribution")
        )

    data = [{"AF": af} for af in af_values]

    chart = (
        alt.Chart(pl.DataFrame(data))
        .mark_bar(color="#556b2f")
        .encode(
            x=alt.X("AF:Q", bin=alt.Bin(maxbins=15), title="AF_max"),
            y=alt.Y("count():Q", title="Count"),
            tooltip=[
                alt.Tooltip("AF:Q", title="AF_max Range", format=".2f"),
                alt.Tooltip("count():Q", title="Count"),
            ],
        )
        .properties(width=200, height=150, title="AF_max Distribution (>0 only)")
    )

    return chart


def create_motif_breakdown_chart(msi_data):
    """Create motif type MSI rate chart"""
    motif_data = []
    for motif_type, data in msi_data["motif_breakdown"].items():
        msi_rate = (data["msi"] / data["total"] * 100) if data["total"] > 0 else 0
        motif_data.append(
            {"Motif": motif_type, "MSI_Rate": msi_rate, "Total": data["total"]}
        )

    chart = (
        alt.Chart(pl.DataFrame(motif_data))
        .mark_bar(color="#8b4513")
        .encode(
            x=alt.X("Motif:N", title="Motif Type", sort="-y"),
            y=alt.Y(
                "MSI_Rate:Q", title="MSI Rate (%)", scale=alt.Scale(domain=[0, 40])
            ),
            tooltip=[
                alt.Tooltip("Motif:N", title="Motif Type"),
                alt.Tooltip("MSI_Rate:Q", title="MSI Rate (%)", format=".1f"),
                alt.Tooltip("Total:Q", title="Total Regions", format=","),
            ],
        )
        .properties(width=200, height=150, title="MSI Rate by Motif Type")
    )

    return chart


def create_msi_score_uncertainty_chart(regional_results):
    """3-bar MSI chart with tooltips - probabilistic, expected, deterministic"""

    try:
        threshold_score = regional_results.get("msi_score", 0.0) or 0.0
        threshold_uncertainty = (
            regional_results.get("msi_score_statistical_uncertainty", 0.0) or 0.0
        )
        deterministic_score = (
            regional_results.get("msi_score_deterministic", 0.0) or 0.0
        )
        expected_score = regional_results.get("msi_score_expected", None)
        expected_uncertainty = regional_results.get(
            "msi_score_expected_uncertainty", None
        )
    except (TypeError, KeyError, AttributeError):
        threshold_score = 0.0
        threshold_uncertainty = 0.0
        deterministic_score = 0.0
        expected_score = None
        expected_uncertainty = None

    chart_data = [
        {
            "Method": "Probabilistic",
            "Score": threshold_score,
            "Uncertainty": (
                f"±{threshold_uncertainty:.3f}%"
                if threshold_uncertainty > 0
                else "±0.000%"
            ),
            "Color": "#8b4513",
        },
        {
            "Method": "Expected",
            "Score": expected_score if expected_score is not None else 0.0,
            "Uncertainty": (
                f"±{expected_uncertainty:.3f}%" if expected_uncertainty else "Pending"
            ),
            "Color": "#daa520" if expected_score is not None else "#cccccc",
        },
        {
            "Method": "Deterministic",
            "Score": deterministic_score,
            "Uncertainty": "±0.000%",
            "Color": "#2f4f4f",
        },
    ]

    try:
        chart = (
            alt.Chart(pl.DataFrame(chart_data))
            .mark_bar(width=50)
            .encode(
                x=alt.X("Method:N", title="Analysis Method"),
                y=alt.Y("Score:Q", title="MSI Score (%)"),
                color=alt.Color(field="Color", type="nominal", scale=None),
                tooltip=[
                    alt.Tooltip("Method:N", title="Method"),
                    alt.Tooltip("Score:Q", title="MSI Score (%)", format=".2f"),
                    alt.Tooltip("Uncertainty:N", title="Uncertainty"),
                ],
            )
            .properties(width=200, height=150, title="MSI Score Comparison")
        )
        return chart
    except Exception:
        return alt.Chart(pl.DataFrame([{"x": 0, "y": 0}])).mark_text(text="Chart Error")


def create_region_count_breakdown_chart(regional_results):
    """
    Region count breakdown
    """

    try:
        prob_unstable = regional_results.get("unstable_regions", 0) or 0
        prob_stable = regional_results.get("probabilistic_stable_regions", 0) or 0

        det_unstable = regional_results.get("regions_with_variants", 0) or 0
        det_stable = regional_results.get("deterministic_stable_regions", 0) or 0

        exp_unstable = regional_results.get("expected_unstable_regions", 0) or 0
        uncertain = regional_results.get("uncertain_regions", 0) or 0

        total_regions = regional_results.get("total_regions", 0) or 0
        exp_stable = max(0, total_regions - exp_unstable - uncertain)

    except (TypeError, KeyError, AttributeError):
        prob_unstable = prob_stable = det_unstable = det_stable = exp_unstable = (
            exp_stable
        ) = uncertain = 0

    chart_data = [
        {
            "Method": "Probabilistic",
            "Type": "Unstable",
            "Count": prob_unstable,
            "Color": "#8b4513",
        },
        {
            "Method": "Probabilistic",
            "Type": "Stable",
            "Count": prob_stable,
            "Color": "#2f4f4f",
        },
        {
            "Method": "Probabilistic",
            "Type": "Uncertain",
            "Count": uncertain,
            "Color": "#daa520",
        },
        {
            "Method": "Expected",
            "Type": "Unstable",
            "Count": exp_unstable,
            "Color": "#8b4513",
        },
        {
            "Method": "Expected",
            "Type": "Stable",
            "Count": exp_stable,
            "Color": "#2f4f4f",
        },
        {
            "Method": "Expected",
            "Type": "Uncertain",
            "Count": uncertain,
            "Color": "#daa520",
        },
        {
            "Method": "Deterministic",
            "Type": "Unstable",
            "Count": det_unstable,
            "Color": "#8b4513",
        },
        {
            "Method": "Deterministic",
            "Type": "Stable",
            "Count": det_stable,
            "Color": "#2f4f4f",
        },
        {
            "Method": "Deterministic",
            "Type": "Uncertain",
            "Count": uncertain,
            "Color": "#daa520",
        },
    ]

    try:
        chart = (
            alt.Chart(pl.DataFrame(chart_data))
            .mark_bar(width=50)
            .encode(
                x=alt.X("Method:N", title="Analysis Method"),
                y=alt.Y("Count:Q", title="Region Count", stack="zero"),
                color=alt.Color(
                    "Type:N", scale=alt.Scale(range=["#8b4513", "#2f4f4f", "#daa520"])
                ),
                tooltip=[
                    alt.Tooltip("Method:N", title="Method"),
                    alt.Tooltip("Type:N", title="Region Type"),
                    alt.Tooltip("Count:Q", title="Count", format=","),
                ],
            )
            .properties(width=200, height=150, title="Region Count Breakdown")
        )
        return chart
    except Exception:
        return alt.Chart(pl.DataFrame([{"x": 0, "y": 0}])).mark_text(text="Chart Error")
