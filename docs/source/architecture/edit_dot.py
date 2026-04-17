import re

# ── Configuration ──────────────────────────────────────────────────────────────

INPUT  = "classes_MMS.dot" # name of the input file
OUTPUT = f"{INPUT[:-4]}_edited.dot" # name of the output file

# Map: substring in node name → fillcolor
COLOR_MAP = {
    # dyn_sys
    "dyn_sys.Forcing"                       : "cyan",
    "dyn_sys.Dynamical_system"              : "lightblue",
    # mms
    "mms.Multiple_scales_system"            : "lightblue",
    "mms.Coord_MMS"                         : "lightgreen",
    "mms.Forcing_MMS"                       : "cyan",
    "mms.Sol_MMS"                           : "pink",
    "mms.Sol_transient"                     : "pink",
    "mms.Substitutions_MMS"                 : "orange",
    # steady_state
    "steady_state.Coord_SS"                 : "lightgreen",
    "steady_state.Forcing_SS"               : "cyan",
    "steady_state.Sol_SS"                   : "pink",
    "steady_state.Sol_bbc"                  : "pink",
    "steady_state.Sol_forced"               : "pink",
    "steady_state.Sol_LC"                   : "pink",
    "steady_state.Stab_SS"                  : "violet",
    "steady_state.Stability"                : "violet",
    "steady_state.Steady_state"             : "lightblue",
    "steady_state.Substitutions_SS"         : "orange",
    # visualisation
    "visualisation.Frequency_response_curve": "lightyellow",
    "visualisation.Amplitude_response_curve": "lightyellow",
    "visualisation.Transient_response"      : "lightyellow",
    "visualisation.Limit_cycle"             : "lightyellow",
}

# Clusters: cluster name → list of substrings identifying member nodes
CLUSTERS = {
    "Core Classes": [
        "dyn_sys.Dynamical_system",
        "mms.Multiple_scales_system",
        "steady_state.Steady_state",
    ],
    "Coordinates Classes": [
        "mms.Coord_MMS",
        "steady_state.Coord_SS",
    ],
    "Solution Classes": [
        "mms.Sol_MMS",
        "mms.Sol_transient",
        "steady_state.Sol_SS",
        "steady_state.Sol_bbc",
        "steady_state.Sol_forced",
        "steady_state.Sol_LC",
    ],
    "Forcing Classes": [
        "dyn_sys.Forcing",
        "mms.Forcing_MMS",
        "steady_state.Forcing_SS",
    ],
    "Substitution Classes": [
        "mms.Substitutions_MMS",
        "steady_state.Substitutions_SS",
    ],
    "Stability Analysis": [
        "steady_state.Stab_SS",
        "steady_state.Stability",
    ],
    "Visualisation Classes": [
        "visualisation.Frequency_response_curve",
        "visualisation.Amplitude_response_curve",
        "visualisation.Transient_response",
        "visualisation.Limit_cycle",             
    ],
}

# ── Read and patch ─────────────────────────────────────────────────────────────

with open(INPUT, "r", encoding="utf-8") as f:
    content = f.read()

# 1. Fix graph name
content = re.sub(r'digraph "[^"]*"', 'digraph "classes_MMS"', content)

# 2. Inject global graph settings after opening brace
global_settings = """rankdir=BT;
splines="curved";
nodesep=0.5;
ranksep=0.75;
charset="utf-8";
fontname="Courier New";
node [fontname="Courier New"];
edge [fontname="Courier New", color="gray", fontcolor="gray", fontsize=10];"""
content = re.sub(r'rankdir=BT\s*\ncharset="utf-8"', global_settings, content)

# 3. Add fillcolor and style="filled" to each node
def patch_node(match):
    node_id = match.group(1)
    attrs   = match.group(2)
    color   = "white"  # default
    for key, val in COLOR_MAP.items():
        if key in node_id:
            color = val
            break
    attrs = attrs.replace('style="solid"', f'style="filled", fillcolor="{color}"')
    return f'"{node_id}" [{attrs}]'

content = re.sub(r'"([^"]+)"\s*\[([^\]]*(?:\[[^\]]*(?:\[[^\]]*\][^\]]*)*\][^\]]*)*)\](?=;)', patch_node, content)

# 4. Patch edge styles
def patch_edge(match):
    src   = match.group(1)
    dst   = match.group(2)
    label = re.search(r'label="([^"]*)"', match.group(3))
    xlabel = f'xlabel="{label.group(1)}"' if label else ''
    return f'"{src}" -> "{dst}" [arrowhead="normal", arrowtail="none", {xlabel}];'

content = re.sub(
    r'"([^"]+)"\s*->\s*"([^"]+)"\s*\[([^\]]+)\];',
    patch_edge,
    content
)

# 5. Inject subgraph clusters before closing brace
clusters_dot = ""
for i, (label, members) in enumerate(CLUSTERS.items()):
    cluster_name = label.lower().replace(" ", "_")
    seen = set()
    node_lines = []
    for node_id in members:
        for n in re.findall(r'"([^"]+)"', content):
            if node_id in n and n not in seen:
                seen.add(n)
                node_lines.append(f'"{n}"')
    nodes = "\n    ".join(node_lines)
    clusters_dot += f"""
  subgraph cluster_{cluster_name} {{
    label="{label}";
    labelloc="b";
    {nodes}
  }}
"""

content = content.rstrip("}\n") + clusters_dot + "\n}\n"

# 6. Truncate long method signatures
def truncate_methods(match):
    label = match.group(0)
    def shorten(m):
        method = m.group(0)
        if len(method) > 50:
            name = method.split('(')[0]
            return name + '(...)<br ALIGN="LEFT"/>'
        return method
    return re.sub(r'\w+\([^<]*\)<br ALIGN="LEFT"/>', shorten, label)

content = re.sub(r'<\{[^}]+\}>', truncate_methods, content)

# ── Write output ───────────────────────────────────────────────────────────────

with open(OUTPUT, "w", encoding="utf-8") as f:
    f.write(content)

print(f"Patched {INPUT} written to {OUTPUT}")
