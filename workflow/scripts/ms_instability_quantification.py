"""
Quantify MSI detection: 
    Ground truth(Mason Variants + Pytrf MSI indels) vs Varlociraptor output
"""

#!/usr/bin/env python3

import argparse
import json
from datetime import datetime

import altair as alt
import pandas as pd
import pysam

alt.data_transformers.disable_max_rows()

def load_ground_truth_variants(vcf_file):
    """
    Load MSI variants from ground truth file (Mason + PyTRF injected variants).
    Returns list of MSI variants with their details.
    """
    msi_variants = []

    try:
        vcf = pysam.VariantFile(vcf_file)

        for record in vcf:
            if record.id and record.id.startswith("MSI_"):
                variant_info = {
                    "chrom": record.chrom,
                    "pos": record.pos,
                    "ref": record.ref,
                    "alt": str(record.alts[0]),
                    "id": record.id,
                    # Extract MSI-specific info
                    "msi_type": record.info.get("MSI_TYPE", "Unknown"),
                    "msi_motif": record.info.get("MSI_MOTIF", "Unknown"),
                    "msi_orig_repeats": record.info.get("MSI_ORIG_REPEATS", 0),
                    "msi_change": record.info.get("MSI_CHANGE", 0),
                }
                msi_variants.append(variant_info)

        vcf.close()

    except FileNotFoundError as e:
        print(f"[MSI-QUANT ERROR] Ground truth file not found: {vcf_file} - {e}")
        raise
    except PermissionError as e:
        print(f"[MSI-QUANT ERROR] Permission denied for ground truth file: {vcf_file} - {e}")
        raise
    except pysam.utils.SamtoolsError as e:
        print(f"[MSI-QUANT ERROR] Invalid VCF format in ground truth file: {e}")
        raise
    except Exception as e:
        print(f"[MSI-QUANT ERROR] Unexpected error reading ground truth: {e}")
        raise

    return msi_variants


def load_detected_variants(vcf_file):
    """
    Load all variants detected by varlociraptor.
    Returns list of all detected variants.
    """
    detected_variants = []

    try:
        vcf = pysam.VariantFile(vcf_file)

        for record in vcf:
            variant_info = {
                "chrom": record.chrom,
                "pos": record.pos,
                "ref": record.ref,
                "alt": str(record.alts[0]),
                "id": record.id,
                # Varlociraptor-specific info
                "svlen": record.info.get("SVLEN", (None,))[0],
                "prob_present": record.info.get("PROB_PRESENT", (None,))[0],
                "prob_artifact": record.info.get("PROB_ARTIFACT", (None,))[0],
            }
            detected_variants.append(variant_info)

        vcf.close()

    except FileNotFoundError:
        print(f"[MSI-QUANT ERROR] Varlociraptor output file not found: {vcf_file}")
        raise
    except PermissionError:
        print(
            f"[MSI-QUANT ERROR] Permission denied for varlociraptor output file: {vcf_file}"
        )
        raise
    except pysam.utils.SamtoolsError as e:
        print(f"[MSI-QUANT ERROR] Invalid VCF format in varlociraptor output: {e}")
        raise
    except Exception as e:
        print(f"[MSI-QUANT ERROR] Unexpected error reading varlociraptor output: {e}")
        raise

    return detected_variants


def compare_msi_variants(ground_truth, detected):
    """
    Compare MSI variants using exact matching.
    """
    detected_map = {(v["chrom"], v["pos"], v["ref"], v["alt"]): v for v in detected}

    matched_msi = []
    missed_msi = []

    for msi_variant in ground_truth:
        key = (
            msi_variant["chrom"],
            msi_variant["pos"],
            msi_variant["ref"],
            msi_variant["alt"],
        )
        if key in detected_map:
            matched_msi.append({**msi_variant, **detected_map[key]})
        else:
            missed_msi.append(msi_variant)

    total_msi = len(ground_truth)
    detected_count = len(matched_msi)
    sensitivity = (detected_count / total_msi * 100) if total_msi > 0 else 0

    total_ins = len([v for v in ground_truth if v["msi_type"] == "INS"])
    total_del = len([v for v in ground_truth if v["msi_type"] == "DEL"])
    detected_ins = len([v for v in matched_msi if v["msi_type"] == "INS"])
    detected_del = len([v for v in matched_msi if v["msi_type"] == "DEL"])

    probs = [v["prob_present"] for v in matched_msi if v["prob_present"] is not None]
    avg_prob = sum(probs) / len(probs) if probs else 0
    high_confidence = len([p for p in probs if p > 0.8])
    low_confidence = len([p for p in probs if p < 0.5])

    return {
        "total_msi": total_msi,
        "detected": detected_count,
        "missed": len(missed_msi),
        "sensitivity": sensitivity,
        "sensitivity_ins": (detected_ins / total_ins * 100) if total_ins > 0 else 0,
        "sensitivity_del": (detected_del / total_del * 100) if total_del > 0 else 0,
        "matched_variants": matched_msi,
        "missed_variants": missed_msi,
        "avg_probability": avg_prob,
        "high_confidence": high_confidence,
        "low_confidence": low_confidence,
    }


def generate_html_report(results, output_path):
    """
    Generate MSI Quantification HTML report.
    """

    matched = results["matched_variants"]
    missed = results["missed_variants"]

    sensitivity_data = [
        {"Type": "Insertions", "Sensitivity": results["sensitivity_ins"]},
        {"Type": "Deletions", "Sensitivity": results["sensitivity_del"]},
    ]

    confidence_data = [
        {"Probability": v["prob_present"]}
        for v in matched
        if v["prob_present"] is not None
    ]

    high_conf_count = len(
        [v for v in matched if v["prob_present"] and v["prob_present"] > 0.8]
    )
    med_conf_count = len(
        [v for v in matched if v["prob_present"] and 0.5 <= v["prob_present"] <= 0.8]
    )
    low_conf_count = len(
        [v for v in matched if v["prob_present"] and v["prob_present"] < 0.5]
    )
    no_conf_count = len([v for v in matched if not v["prob_present"]])

    confidence_breakdown = [
        {"Level": "High (>0.8)", "Count": high_conf_count},
        {"Level": "Medium (0.5-0.8)", "Count": med_conf_count},
        {"Level": "Low (<0.5)", "Count": low_conf_count},
        {"Level": "No Data", "Count": no_conf_count},
    ]

    sensitivity_chart = (
        alt.Chart(pd.DataFrame(sensitivity_data))
        .mark_bar(color="#8b4513", width=50)
        .encode(
            x=alt.X("Type:N", title="MSI Type"),
            y=alt.Y(
                "Sensitivity:Q",
                scale=alt.Scale(domain=[0, 100]),
                title="Sensitivity (%)",
            ),
            tooltip=[
                alt.Tooltip("Type:N", title="Type"),
                alt.Tooltip("Sensitivity:Q", title="Sensitivity (%)", format=".1f"),
            ],
        )
        .properties(width=200, height=150, title="Sensitivity by Type")
    )

    if confidence_data:
        confidence_chart = (
            alt.Chart(pd.DataFrame(confidence_data))
            .mark_circle(size=100, color="#2f4f4f", opacity=0.8)
            .encode(
                x=alt.X("Probability:Q", title="Detection Probability"),
                y=alt.value(75),
                tooltip=[
                    alt.Tooltip("Probability:Q", title="Probability", format=".3f")
                ],
            )
            .properties(width=200, height=150, title="Confidence Values")
        )
    else:
        confidence_chart = (
            alt.Chart(pd.DataFrame([{"x": 0, "y": 0}]))
            .mark_text(text="No data", size=14)
            .encode(x="x:Q", y="y:Q")
            .properties(width=200, height=150, title="Confidence Values")
        )

    conf_breakdown_chart = (
        alt.Chart(pd.DataFrame(confidence_breakdown))
        .mark_bar(color="#556b2f", width=30)
        .encode(
            x=alt.X("Level:N", title="Confidence Level"),
            y=alt.Y("Count:Q", title="Number of Variants"),
            tooltip=[
                alt.Tooltip("Level:N", title="Confidence Level"),
                alt.Tooltip("Count:Q", title="Count"),
            ],
        )
        .properties(width=200, height=150, title="Detection by Confidence")
    )

    table_data = []
    for v in matched:
        conf_level = (
            "high"
            if v["prob_present"] and v["prob_present"] > 0.8
            else (
                "medium"
                if v["prob_present"] and 0.5 <= v["prob_present"] <= 0.8
                else "low" if v["prob_present"] and v["prob_present"] < 0.5 else "none"
            )
        )

        table_data.append(
            {
                "pos": f"{v['chrom']}:{v['pos']}",
                "change": f"{v['ref']}→{v['alt']}",
                "type": v["msi_type"],
                "motif": v["msi_motif"],
                "prob": (
                    f"{v['prob_present']:.3f}" if v["prob_present"] is not None else "—"
                ),
                "status": "Detected",
                "confidence": conf_level,
            }
        )

    for v in missed:
        table_data.append(
            {
                "pos": f"{v['chrom']}:{v['pos']}",
                "change": f"{v['ref']}→{v['alt']}",
                "type": v["msi_type"],
                "motif": v["msi_motif"],
                "prob": "—",
                "status": "Missed",
                "confidence": "none",
            }
        )

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
            border-bottom: 1px solid #666;
        }}
        .header-info {{
            text-align: center;
            font-style: italic;
            margin-bottom: 40px;
            color: #666;
        }}
        .metrics {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
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
            grid-template-columns: 1fr 1fr 1fr;
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
        .status-detected {{
            color: #2f4f4f;
            font-weight: bold;
        }}
        .status-missed {{
            color: #8b4513;
            font-weight: bold;
        }}
        .filter-btn {{
            padding: 8px 16px;
            margin: 5px;
            border: 1px solid #666;
            background: #fefdfb;
            cursor: pointer;
            font-family: inherit;
        }}
        .filter-btn.active {{
            background: #333;
            color: white;
        }}
        .pagination {{
            display: flex;
            justify-content: center;
            align-items: center;
            margin: 20px 0;
            gap: 10px;
        }}
        .pagination button {{
            padding: 8px 12px;
            border: 1px solid #666;
            background: #fefdfb;
            cursor: pointer;
        }}
        .pagination button:disabled {{
            opacity: 0.5;
            cursor: not-allowed;
        }}
        .pagination button.active {{
            background: #333;
            color: white;
        }}
    </style>
</head>
<body>
    <h1>MSI Quantification Report</h1>

    <div class="header-info">
        Microsatellite Instability Detection Analysis<br>
        Generated {timestamp}
    </div>

    <div class="metrics">
        <div class="metric">
            <span class="metric-value">{results['sensitivity']:.1f}%</span>
            <div class="metric-label">Overall Sensitivity</div>
        </div>
        <div class="metric">
            <span class="metric-value">{results['detected']}</span>
            <div class="metric-label">MSI Variants Detected</div>
        </div>
        <div class="metric">
            <span class="metric-value">{results['missed']}</span>
            <div class="metric-label">MSI Variants Missed</div>
        </div>
        <div class="metric">
            <span class="metric-value">{results['total_msi']}</span>
            <div class="metric-label">Ground Truth Total</div>
        </div>
        <div class="metric">
            <span class="metric-value">{results['avg_probability']:.2f}</span>
            <div class="metric-label">Average Confidence</div>
        </div>
        <div class="metric">
            <span class="metric-value">{high_conf_count}</span>
            <div class="metric-label">High Confidence (>0.8)</div>
        </div>
    </div>

    <div class="charts">
        <div class="chart-container">
            <div id="sensitivity-chart"></div>
        </div>
        <div class="chart-container">
            <div id="confidence-chart"></div>
        </div>
        <div class="chart-container">
            <div id="breakdown-chart"></div>
        </div>
    </div>

    <h2>MSI Variant Details</h2>
    
    <div style="margin-bottom: 15px;">
        <button class="filter-btn active" onclick="filterTable('all')" id="btn-all">
            All ({len(table_data)})
        </button>
        <button class="filter-btn" onclick="filterTable('detected')" id="btn-detected">
            Detected ({results['detected']})
        </button>
        <button class="filter-btn" onclick="filterTable('missed')" id="btn-missed">
            Missed ({results['missed']})
        </button>
        <button class="filter-btn" onclick="filterTable('high')" id="btn-high">
            High Conf ({high_conf_count})
        </button>
        <button class="filter-btn" onclick="filterTable('medium')" id="btn-medium">
            Med Conf ({med_conf_count})
        </button>
        <button class="filter-btn" onclick="filterTable('low')" id="btn-low">
            Low Conf ({low_conf_count})
        </button>
    </div>

    <table>
        <thead>
            <tr>
                <th>Position</th>
                <th>Sequence Change</th>
                <th>Type</th>
                <th>Motif</th>
                <th>Probability</th>
                <th>Status</th>
            </tr>
        </thead>
        <tbody id="variant-table">
        </tbody>
    </table>

    <div class="pagination">
        <button onclick="changePage(-1)" id="prev-btn">← Previous</button>
        <span id="page-info">Page 1 of 1</span>
        <button onclick="changePage(1)" id="next-btn">Next →</button>
        <select id="page-size" onchange="changePageSize()">
            <option value="20">20 per page</option>
            <option value="50">50 per page</option>
            <option value="100">100 per page</option>
        </select>
    </div>

    <script>
        // Chart data
        const sensitivityChart = {sensitivity_chart.to_json()};
        const confidenceChart = {confidence_chart.to_json()};
        const breakdownChart = {conf_breakdown_chart.to_json()};
        
        // Render charts
        vegaEmbed('#sensitivity-chart', sensitivityChart, {{actions: false}});
        vegaEmbed('#confidence-chart', confidenceChart, {{actions: false}});
        vegaEmbed('#breakdown-chart', breakdownChart, {{actions: false}});
        
        // Table data and pagination
        const variants = {json.dumps(table_data)};
        let currentFilter = 'all';
        let currentPage = 1;
        let pageSize = 20;
        
        function getFilteredVariants() {{
            if (currentFilter === 'all') return variants;
            if (currentFilter === 'detected') return variants.filter(v => v.status === 'Detected');
            if (currentFilter === 'missed') return variants.filter(v => v.status === 'Missed');
            if (currentFilter === 'high') return variants.filter(v => v.confidence === 'high');
            if (currentFilter === 'medium') return variants.filter(v => v.confidence === 'medium');
            if (currentFilter === 'low') return variants.filter(v => v.confidence === 'low');
            return variants;
        }}

        function renderTable() {{
            const filtered = getFilteredVariants();
            const totalPages = Math.ceil(filtered.length / pageSize);
            const startIdx = (currentPage - 1) * pageSize;
            const endIdx = startIdx + pageSize;
            const pageData = filtered.slice(startIdx, endIdx);

            const tbody = document.getElementById('variant-table');
            tbody.innerHTML = pageData.map(v => `
                <tr>
                    <td>${{v.pos}}</td>
                    <td>${{v.change}}</td>
                    <td>${{v.type}}</td>
                    <td>${{v.motif}}</td>
                    <td>${{v.prob}}</td>
                    <td class="status-${{v.status.toLowerCase()}}">${{v.status}}</td>
                </tr>
            `).join('');
            
            // Update pagination
            document.getElementById('page-info').textContent = 
                `Page ${{currentPage}} of ${{totalPages}} (${{filtered.length}} variants)`;
            document.getElementById('prev-btn').disabled = currentPage === 1;
            document.getElementById('next-btn').disabled = currentPage === totalPages;
        }}

        function filterTable(filter) {{
            currentFilter = filter;
            currentPage = 1;
            document.querySelectorAll('.filter-btn').forEach(btn => btn.classList.remove('active'));
            document.getElementById(`btn-${{filter}}`).classList.add('active');
            renderTable();
        }}

        function changePage(direction) {{
            const filtered = getFilteredVariants();
            const totalPages = Math.ceil(filtered.length / pageSize);
            const newPage = currentPage + direction;
            if (newPage >= 1 && newPage <= totalPages) {{
                currentPage = newPage;
                renderTable();
            }}
        }}

        function changePageSize() {{
            pageSize = parseInt(document.getElementById('page-size').value);
            currentPage = 1;
            renderTable();
        }}

        renderTable();
    </script>
</body>
</html>
    """

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"HTML report saved: {output_path}")


def main():
    """
    Main function to parse arguments and run the MSI quantification.
    """
    parser = argparse.ArgumentParser(description="Quantify MSI detection")
    parser.add_argument("--ground-truth", required=True, help="Injected MSI variants")
    parser.add_argument(
        "--varlociraptor-output", required=True, help="Varlociraptor output"
    )
    parser.add_argument("--html-output", help="Path to save HTML results page")

    args = parser.parse_args()
    print("Starting MSI detection analysis...")

    # Load variants from both files
    print("\n=== LOADING DATA ===" + "=" * 60)
    ground_truth = load_ground_truth_variants(args.ground_truth)
    detected = load_detected_variants(args.varlociraptor_output)
    print(f"Loaded {len(ground_truth)} MSI variants from ground truth")
    print(f"Loaded {len(detected)} total variants from varlociraptor")

    results = compare_msi_variants(ground_truth, detected)
    print(
        f"Sensitivity: {results['sensitivity']:.1f}% ({results['detected']}/{results['total_msi']})"
    )
    print("=" * 80 + "\n")

    if args.html_output:
        generate_html_report(results, args.html_output)


if __name__ == "__main__":
    main()
