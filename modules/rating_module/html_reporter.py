#!/usr/bin/env python3
import os
import json
import argparse

def generate_html(eval_json_path, merge_json_path, out_html_path):
    # 1. Load Data Safely
    eval_data = {}
    if os.path.exists(eval_json_path):
        with open(eval_json_path, 'r') as f:
            eval_data = json.load(f)
            
    merge_data = {}
    if os.path.exists(merge_json_path):
        with open(merge_json_path, 'r') as f:
            merge_data = json.load(f)

    if not eval_data and not merge_data:
        print("[Error] Both JSON files are missing or empty.")
        return

    # 2. Extract Global Stats
    total_seqs = merge_data.get("Total_Sequences", 0)
    tools = list(eval_data.keys())

    # Helper function to extract specific TE stats safely
    def get_te_stat(tool_data, te_class, metric):
        try:
            if metric == 'count':
                return tool_data["4_Post_Repair_Yield"]["Distribution"].get(te_class, {}).get("count", 0)
            elif metric in ['Median', 'N50']:
                return tool_data["5_Post_Repair_Lengths"].get(te_class, {}).get(metric, 0)
        except:
            return 0
        return 0

    # 3. Build HTML Structure
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en" class="light">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>RepLibBuilder v1.0.1 Report</title>
        <script src="https://cdn.tailwindcss.com"></script>
        <script>
            tailwind.config = {{
                darkMode: 'class',
                theme: {{
                    extend: {{}}
                }}
            }}
        </script>
        <script src="https://cdn.jsdelivr.net/npm/echarts@5.5.0/dist/echarts.min.js"></script>
        <style>
            body {{ transition: background-color 0.3s, color 0.3s; }}
            .card {{ border-radius: 0.5rem; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1); padding: 1.5rem; transition: background-color 0.3s, border-color 0.3s; }}
            .table-row-even {{ background-color: #f9fafb; }}
            .dark .table-row-even {{ background-color: #1f2937; }}
        </style>
    </head>
    <body class="bg-gray-50 text-gray-800 dark:bg-gray-900 dark:text-gray-200 font-sans antialiased p-6">

        <div class="max-w-7xl mx-auto">
            <header class="mb-8 border-b border-gray-200 dark:border-gray-700 pb-6 flex justify-between items-end">
                <div>
                    <h1 class="text-4xl font-extrabold text-gray-900 dark:text-white">RepLibBuilder v1.0.1 Evaluation Report</h1>
                    <p class="text-gray-500 dark:text-gray-400 mt-2 text-lg">Automated diagnostics, standardization, and quality assessment of Repeat Libraries.</p>
                </div>
                <button id="themeToggle" class="px-4 py-2 bg-indigo-100 text-indigo-800 dark:bg-indigo-900 dark:text-indigo-200 rounded-full font-semibold shadow-sm hover:shadow transition-shadow">
                    🌙 Dark Mode
                </button>
            </header>

            <div class="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-8">
                <div class="card bg-blue-50 border border-blue-100 dark:bg-blue-900/20 dark:border-blue-800/50">
                    <h2 class="text-lg font-bold text-blue-800 dark:text-blue-300 mb-2">📊 Confidence Index Scoring</h2>
                    <div class="text-sm text-blue-900 dark:text-blue-200 space-y-2">
                        <ul class="list-disc pl-5 space-y-1">
                            <li><strong>Track A (Domain Evidence):</strong> Perfect &mdash; both sides confident &amp; identical at sf (+3), Asymmetric &mdash; at least one side sf=Unknown, no refutation possible (+2), Fuzzy &mdash; sf-level disagreement or unilateral rescue (+1), Missed detection (-1), Major class conflict (-2).</li>
                            <li><strong>Track B (Structural Heuristics):</strong> Valid short SINE/MITEs (+0.5), Truncated elements (0), Unknown/Unverified (0).</li>
                        </ul>
                        <div class="mt-3 p-3 bg-white/60 dark:bg-black/20 rounded shadow-sm space-y-2 border border-blue-200 dark:border-blue-700/50">
                            <p><strong>Absolute Score:</strong> The raw cumulative score of the entire library.</p>
                            <p><strong>Normalized Score:</strong> (Absolute Score / Total Sequences) &times; 100. Allows size-independent comparison across different tools.</p>
                            <p><strong>Mean Scored:</strong> Absolute Score divided by the number of sequences with a non-zero score.</p>
                        </div>
                    </div>
                </div>

                <div class="card bg-green-50 border border-green-100 dark:bg-green-900/20 dark:border-green-800/50">
                    <h2 class="text-lg font-bold text-green-800 dark:text-green-300 mb-2">🛠️ Active Repair Tags Definition</h2>
                    <div class="text-sm text-green-900 dark:text-green-200 grid grid-cols-1 sm:grid-cols-2 gap-x-4">
                        <ul class="list-disc pl-5 space-y-1">
                            <li><strong>Confirmed:</strong> Domain evidence perfectly matches prediction or completes a subclass.</li>
                            <li><strong>Corrected:</strong> Major conflict overridden by strong domain evidence.</li>
                            <li><strong>Recovered:</strong> Originally "Unknown", successfully identified via domain.</li>
                        </ul>
                        <ul class="list-disc pl-5 space-y-1">
                            <li><strong>Exempted:</strong> Valid short structural TEs lacking domains (e.g., SINE < 1kb).</li>
                            <li><strong>Unverified:</strong> No domain detected; retained as-is. Covers both normal-length elements lacking a domain and short (&lt; 1 kb) truncated copies of autonomous classes (LTR/LINE/DNA/RC/PLE), whose absent domain is expected and not penalized.</li>
                        </ul>
                    </div>
                </div>
            </div>

            <h2 class="text-2xl font-semibold mb-4 text-gray-800 dark:text-gray-100 border-l-4 border-indigo-500 pl-3">1. Global Pool Summary</h2>
            <div class="grid grid-cols-1 md:grid-cols-3 gap-6 mb-10">
                <div class="card bg-gradient-to-br from-indigo-500 to-purple-600 border-0 text-white flex flex-col justify-center items-center shadow-lg">
                    <h3 class="text-lg font-medium opacity-90">Total Pooled Sequences</h3>
                    <p class="text-6xl font-black mt-2 drop-shadow-md">{total_seqs:,}</p>
                </div>
                <div class="card bg-white dark:bg-gray-800 border dark:border-gray-700 md:col-span-2">
                    <h3 class="text-lg font-medium text-gray-600 dark:text-gray-300 mb-2">Sequence Source Composition</h3>
                    <div id="globalPieChart" style="width: 100%; height: 250px;"></div>
                </div>
            </div>
    """

    # Section 2: De Novo Diagnostics (Only if tools exist)
    if tools:
        html_content += f"""
            <h2 class="text-2xl font-semibold mb-4 text-gray-800 dark:text-gray-100 border-l-4 border-indigo-500 pl-3">2. Active Repair & Tag Distribution</h2>
            <div class="card bg-white dark:bg-gray-800 border dark:border-gray-700 mb-10">
                <h3 class="text-md font-medium text-gray-500 dark:text-gray-400 mb-4">Cross-Software Diagnostic Tags</h3>
                <div id="tagsBarChart" style="width: 100%; height: 350px;"></div>
            </div>
            
            <h2 class="text-2xl font-semibold mb-4 text-gray-800 dark:text-gray-100 border-l-4 border-indigo-500 pl-3">3. Tool-Specific Quality Profiling</h2>
        """
        
        target_classes = ['LTR', 'LINE', 'SINE', 'DNA', 'Unknown']

        for tool in tools:
            ci = eval_data[tool]["1_Confidence_Index"]
            abs_score = ci["Absolute_Score"]
            norm_score = ci["Normalized_Score"]
            mean_score = ci["Mean_Scored"]
            
            # Generate Table Rows for specific TE Classes
            table_rows = ""
            for i, te_cls in enumerate(target_classes):
                count = get_te_stat(eval_data[tool], te_cls, 'count')
                median = get_te_stat(eval_data[tool], te_cls, 'Median')
                n50 = get_te_stat(eval_data[tool], te_cls, 'N50')
                row_class = "table-row-even" if i % 2 == 0 else "bg-white dark:bg-gray-800"
                
                table_rows += f"""
                <tr class="{row_class} border-b dark:border-gray-700 text-center text-sm">
                    <td class="py-3 px-4 font-semibold text-gray-700 dark:text-gray-300">{te_cls}</td>
                    <td class="py-3 px-4 text-gray-600 dark:text-gray-400">{count:,}</td>
                    <td class="py-3 px-4 text-gray-600 dark:text-gray-400">{median:,.1f}</td>
                    <td class="py-3 px-4 font-medium text-indigo-600 dark:text-indigo-400">{n50:,}</td>
                </tr>
                """

            html_content += f"""
            <div class="mb-12 bg-white dark:bg-gray-800 rounded-lg shadow-sm border border-gray-200 dark:border-gray-700 overflow-hidden">
                <div class="bg-gray-50 dark:bg-gray-900/50 px-6 py-4 border-b border-gray-200 dark:border-gray-700">
                    <h3 class="text-xl font-bold text-gray-800 dark:text-gray-100">{tool} Performance Profile</h3>
                </div>
                
                <div class="p-6">
                    <div class="grid grid-cols-1 md:grid-cols-3 gap-6 mb-8">
                        <div class="bg-blue-50 dark:bg-blue-900/20 border-l-4 border-blue-500 p-4 rounded-r-lg">
                            <p class="text-sm text-blue-600 dark:text-blue-400 font-bold uppercase tracking-wider">Normalized Score</p>
                            <p class="text-4xl font-black text-blue-900 dark:text-blue-100 mt-1">{norm_score}</p>
                            <p class="text-xs text-blue-500 dark:text-blue-300 mt-1">Per 100 sequences</p>
                        </div>
                        <div class="bg-indigo-50 dark:bg-indigo-900/20 border-l-4 border-indigo-500 p-4 rounded-r-lg">
                            <p class="text-sm text-indigo-600 dark:text-indigo-400 font-bold uppercase tracking-wider">Mean Scored</p>
                            <p class="text-4xl font-black text-indigo-900 dark:text-indigo-100 mt-1">{mean_score}</p>
                            <p class="text-xs text-indigo-500 dark:text-indigo-300 mt-1">Score / non-zero seqs</p>
                        </div>
                        <div class="bg-purple-50 dark:bg-purple-900/20 border-l-4 border-purple-500 p-4 rounded-r-lg">
                            <p class="text-sm text-purple-600 dark:text-purple-400 font-bold uppercase tracking-wider">Absolute Score</p>
                            <p class="text-4xl font-black text-purple-900 dark:text-purple-100 mt-1">{abs_score}</p>
                            <p class="text-xs text-purple-500 dark:text-purple-300 mt-1">Total raw points</p>
                        </div>
                    </div>

                    <div class="grid grid-cols-1 lg:grid-cols-2 gap-8">
                        <div>
                            <h4 class="text-sm font-bold text-gray-500 dark:text-gray-400 uppercase tracking-wider mb-3">Key Structural Elements</h4>
                            <div class="overflow-x-auto border border-gray-200 dark:border-gray-700 rounded-lg">
                                <table class="min-w-full">
                                    <thead class="bg-gray-100 dark:bg-gray-900 border-b border-gray-200 dark:border-gray-700 text-gray-600 dark:text-gray-300 text-sm">
                                        <tr>
                                            <th class="py-3 px-4 text-center font-semibold">TE Class</th>
                                            <th class="py-3 px-4 text-center font-semibold">Count</th>
                                            <th class="py-3 px-4 text-center font-semibold">Median (bp)</th>
                                            <th class="py-3 px-4 text-center font-semibold border-l border-gray-200 dark:border-gray-700 text-indigo-700 dark:text-indigo-400">N50 (bp)</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        {table_rows}
                                    </tbody>
                                </table>
                            </div>
                        </div>

                        <div>
                            <h4 class="text-sm font-bold text-gray-500 dark:text-gray-400 uppercase tracking-wider mb-3 text-center">Diagnostic Signatures (All Rules)</h4>
                            <div id="radarChart_{tool}" style="width: 100%; height: 350px;"></div>
                        </div>
                    </div>
                </div>
            </div>
            """

    # 4. Inject JavaScript for ECharts and Dark Mode Toggle
    html_content += f"""
        </div> <script>
            const mergeData = {json.dumps(merge_data)};
            const evalData = {json.dumps(eval_data)};
            const tools = {json.dumps(tools)};

            let pieChart, tagsChart;
            const radarCharts = {{}};

            function initCharts(isDark) {{
                const textColor = isDark ? '#e5e7eb' : '#374151';
                
                // 1. Global Pie Chart
                if (document.getElementById('globalPieChart')) {{
                    if(pieChart) pieChart.dispose();
                    pieChart = echarts.init(document.getElementById('globalPieChart'));
                    let pieData = [];
                    if (mergeData.Public_DB) {{
                        for (const [key, val] of Object.entries(mergeData.Public_DB)) {{
                            if (val > 0) pieData.push({{name: key + ' (Public)', value: val}});
                        }}
                    }}
                    if (mergeData.Denovo_DB) {{
                        for (const [key, val] of Object.entries(mergeData.Denovo_DB)) {{
                            let sum = Object.values(val).reduce((a, b) => a + b, 0);
                            if (sum > 0) pieData.push({{name: key + ' (Denovo)', value: sum}});
                        }}
                    }}
                    
                    pieChart.setOption({{
                        textStyle: {{ color: textColor }},
                        tooltip: {{ trigger: 'item' }},
                        legend: {{ orient: 'vertical', left: 'left', textStyle: {{ color: textColor }} }},
                        series: [{{
                            type: 'pie',
                            radius: ['40%', '75%'],
                            itemStyle: {{ borderRadius: 8, borderColor: isDark ? '#1f2937' : '#fff', borderWidth: 2 }},
                            label: {{ color: textColor }},
                            data: pieData
                        }}]
                    }});
                }}

                // 2. Tags Stacked Bar Chart
                if (document.getElementById('tagsBarChart') && tools.length > 0) {{
                    if(tagsChart) tagsChart.dispose();
                    tagsChart = echarts.init(document.getElementById('tagsBarChart'));
                    const tags = ['Confirmed', 'Corrected', 'Recovered', 'Exempted', 'Unverified'];
                    const seriesData = tags.map(tag => ({{
                        name: tag,
                        type: 'bar',
                        stack: 'total',
                        emphasis: {{ focus: 'series' }},
                        data: tools.map(tool => mergeData.Denovo_DB[tool][tag] || 0)
                    }}));

                    tagsChart.setOption({{
                        textStyle: {{ color: textColor }},
                        tooltip: {{ trigger: 'axis', axisPointer: {{ type: 'shadow' }} }},
                        legend: {{ data: tags, bottom: 0, textStyle: {{ color: textColor }} }},
                        grid: {{ left: '3%', right: '4%', bottom: '10%', containLabel: true }},
                        xAxis: {{ type: 'value', splitLine: {{ lineStyle: {{ color: isDark ? '#374151' : '#e5e7eb' }} }} }},
                        yAxis: {{ type: 'category', data: tools }},
                        color: ['#10b981', '#3b82f6', '#8b5cf6', '#f59e0b', '#cbd5e1'],
                        series: seriesData
                    }});
                }}

                // 3. Radar Charts Setup (All 8 Rules)
                let radarMax = {{ perfect: 10, asymmetric: 10, fuzzy: 10, miss: 10, conflict: 10, sineExempt: 10, truncExempt: 10, unknown: 10 }};
                if (tools.length > 0) {{
                    tools.forEach(t => {{
                        let st = evalData[t]["1_Confidence_Index"]["Detailed_Stats"];
                        if ((st.TrackA_Perfect_Plus3 || 0) > radarMax.perfect) radarMax.perfect = st.TrackA_Perfect_Plus3;
                        if ((st.TrackA_Asymmetric_Plus2 || 0) > radarMax.asymmetric) radarMax.asymmetric = st.TrackA_Asymmetric_Plus2;
                        if ((st.TrackA_Fuzzy_Plus1 || 0) > radarMax.fuzzy) radarMax.fuzzy = st.TrackA_Fuzzy_Plus1;
                        if ((st.TrackA_Miss_Minus1 || 0) > radarMax.miss) radarMax.miss = st.TrackA_Miss_Minus1;
                        if ((st.TrackA_Conflict_Minus2 || 0) > radarMax.conflict) radarMax.conflict = st.TrackA_Conflict_Minus2;
                        if ((st.TrackB_SINE_Exempt_Plus0_5 || 0) > radarMax.sineExempt) radarMax.sineExempt = st.TrackB_SINE_Exempt_Plus0_5;
                        if ((st.TrackB_Trunc_Exempt_0 || 0) > radarMax.truncExempt) radarMax.truncExempt = st.TrackB_Trunc_Exempt_0;
                        if ((st.TrackB_Unknown_0 || 0) > radarMax.unknown) radarMax.unknown = st.TrackB_Unknown_0;
                    }});
                    Object.keys(radarMax).forEach(k => radarMax[k] = Math.ceil(radarMax[k] * 1.1));
                }}

                tools.forEach(tool => {{
                    if (document.getElementById('radarChart_' + tool)) {{
                        if(radarCharts[tool]) radarCharts[tool].dispose();
                        const rc = echarts.init(document.getElementById('radarChart_' + tool));
                        radarCharts[tool] = rc;
                        
                        const stats = evalData[tool]["1_Confidence_Index"]["Detailed_Stats"];
                        
                        rc.setOption({{
                            tooltip: {{ trigger: 'item' }},
                            radar: {{
                                indicator: [
                                    {{ name: 'Perfect (+3)', max: radarMax.perfect }},
                                    {{ name: 'Asymmetric (+2)', max: radarMax.asymmetric }},
                                    {{ name: 'Fuzzy (+1)', max: radarMax.fuzzy }},
                                    {{ name: 'Miss (-1)', max: radarMax.miss }},
                                    {{ name: 'Conflict (-2)', max: radarMax.conflict }},
                                    {{ name: 'SINE/MITE Exmpt (+0.5)', max: radarMax.sineExempt }},
                                    {{ name: 'Trunc Exmpt (0)', max: radarMax.truncExempt }},
                                    {{ name: 'Unknown (0)', max: radarMax.unknown }}
                                ],
                                radius: '60%', 
                                splitNumber: 4,
                                splitLine: {{ lineStyle: {{ color: isDark ? '#374151' : '#e5e7eb' }} }},
                                splitArea: {{ show: false }},
                                axisName: {{ color: isDark ? '#9ca3af' : '#64748b', fontSize: 10, padding: [3, 5] }}
                            }},
                            series: [{{
                                type: 'radar',
                                areaStyle: {{ color: 'rgba(99, 102, 241, 0.2)' }},
                                lineStyle: {{ color: '#4f46e5', width: 2 }},
                                itemStyle: {{ color: '#4f46e5' }},
                                data: [{{
                                    value: [
                                        stats.TrackA_Perfect_Plus3 || 0,
                                        stats.TrackA_Asymmetric_Plus2 || 0,
                                        stats.TrackA_Fuzzy_Plus1 || 0,
                                        stats.TrackA_Miss_Minus1 || 0,
                                        stats.TrackA_Conflict_Minus2 || 0,
                                        stats.TrackB_SINE_Exempt_Plus0_5 || 0,
                                        stats.TrackB_Trunc_Exempt_0 || 0,
                                        stats.TrackB_Unknown_0 || 0
                                    ],
                                    name: tool
                                }}]
                            }}]
                        }});
                    }}
                }});
            }}

            // Theme Toggle Logic
            const htmlEl = document.documentElement;
            const themeBtn = document.getElementById('themeToggle');
            
            // Check system preference
            if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {{
                htmlEl.classList.add('dark');
            }}

            initCharts(htmlEl.classList.contains('dark'));

            themeBtn.addEventListener('click', () => {{
                htmlEl.classList.toggle('dark');
                const isDark = htmlEl.classList.contains('dark');
                themeBtn.innerHTML = isDark ? '☀️ Light Mode' : '🌙 Dark Mode';
                initCharts(isDark); // Re-init charts to update colors smoothly
            }});

            // Resize handler
            window.addEventListener('resize', function() {{
                if(pieChart) pieChart.resize();
                if(tagsChart) tagsChart.resize();
                tools.forEach(t => {{ if(radarCharts[t]) radarCharts[t].resize(); }});
            }});
        </script>
    </body>
    </html>
    """

    with open(out_html_path, 'w') as f:
        f.write(html_content)
    
    print(f"\n[Success] Interactive HTML report generated at: {out_html_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standalone visualizer for RepLibBuilder JSON reports.")
    parser.add_argument("--eval", required=True, help="Path to evaluation_report.json")
    parser.add_argument("--merge", required=True, help="Path to merge_stats.json")
    parser.add_argument("--out", default="evaluation_report.html", help="Output HTML path")
    
    args = parser.parse_args()
    generate_html(args.eval, args.merge, args.out)
