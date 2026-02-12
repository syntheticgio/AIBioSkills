# CLAUDE.md — Project Instructions

## What This Project Is

An MCP server for bioinformatics (`datasets-mcp-server`) exposing tools backed by NCBI Datasets CLI, BioPython, Ensembl REST API, and UniProt REST API.

## Build & Test

```bash
pip install -e .                                           # install
python -c "from datasets_mcp.server import app"           # verify imports
npx @modelcontextprotocol/inspector datasets-mcp-server   # interactive test
```

## Code Layout

| Path | Purpose |
|------|---------|
| `src/datasets_mcp/server.py` | All tool definitions (`list_tools`) and handlers (`call_tool`) |
| `pyproject.toml` | Dependencies and entry point |
| `.claude/commands/project-summary.md` | `/project-summary` skill definition |
| `.claude/commands/gene-report.md` | `/gene-report` skill — multi-database gene report |
| `.claude/commands/variant-report.md` | `/variant-report` skill — variant annotation report |
| `.claude/commands/protein-report.md` | `/protein-report` skill — protein-centric report |
| `.claude/commands/pathway-report.md` | `/pathway-report` skill — pathway deep-dive report |
| `.claude/commands/drug-target-report.md` | `/drug-target-report` skill — druggability assessment |
| `.claude/commands/literature-agent.md` | `/literature-agent` agent — PubMed literature search & summary |
| `.claude/commands/comparative-genomics-agent.md` | `/comparative-genomics-agent` agent — multi-species gene comparison |
| `.claude/commands/clinical-variant-agent.md` | `/clinical-variant-agent` agent — full clinical variant workup |
| `.claude/commands/gene-list-agent.md` | `/gene-list-agent` agent — gene list functional analysis |
| `.claude/commands/structure-agent.md` | `/structure-agent` agent — structural biology deep-dive |
| `.claude/commands/statistical-methods-agent.md` | `/statistical-methods-agent` agent — paper statistical analysis |

## Documentation Files — Keep In Sync

When you add, remove, or modify a tool, skill, or agent, you **must** update every applicable file below. Do not treat any single file as the source of truth — they serve different audiences and purposes.

### File Inventory

| File | Purpose | What to update |
|------|---------|----------------|
| **`README.md`** | User-facing overview, setup, and usage examples | Add the tool to the Features list (correct section). Add 1-2 natural-language example prompts and a JSON parameter block under Usage Examples. Update Planned Skills/Agents if promoting from planned to implemented. |
| **`EXAMPLES.md`** | Extended usage examples and workflow patterns | Add example queries showing how a user would invoke the tool via natural language. Include the JSON parameters the assistant would send. Add workflow examples if the tool combines well with others. |
| **`PROJECT_INVENTORY.md`** | Technical reference — every tool with its full parameter table | Add a subsection under the correct category header. Include a parameter table (Parameter / Type / Required / Description). Document the return value. |
| **`ROADMAP.md`** | Tracking document for implemented and planned items | Move the tool from "Planned" to "Implemented" (or add it to "Implemented" if it was never planned). Update planned skills/agents if affected. |
| **`.claude/commands/project-summary.md`** | `/project-summary` skill prompt | Update the example list and tool groupings if the new tool changes what a summary should show. |

### Checklist for Adding a New Tool

1. **`server.py`** — Add the tool definition in `list_tools()` and its handler in `call_tool()`. Place it in the correct section (NCBI / BioPython / Ensembl / UniProt / new section). If the tool doesn't need the NCBI CLI, add its name to `REST_API_TOOLS` or `BIOPYTHON_TOOLS`.
2. **`pyproject.toml`** — Add any new dependencies.
3. **`PROJECT_INVENTORY.md`** — Add the full parameter table and return description.
4. **`README.md`** — Add to Features list and Usage Examples (with natural-language prompts + JSON).
5. **`EXAMPLES.md`** — Add example queries, including multi-tool workflow examples if relevant.
6. **`ROADMAP.md`** — Add to "Implemented Tools" table. Remove from "Planned" if applicable.
7. **`.claude/commands/project-summary.md`** — Update if the new tool changes the categories or example set.

### Checklist for Adding a New Skill

1. **`.claude/commands/<skill-name>.md`** — Create the skill definition file.
2. **`PROJECT_INVENTORY.md`** — Add under "Custom Skills" with a description.
3. **`README.md`** — Move from "Planned Skills" to a new "Skills" section (or add directly). Include an example prompt.
4. **`EXAMPLES.md`** — Add usage examples showing the slash command.
5. **`ROADMAP.md`** — Move from "Planned Skills" to a new "Implemented Skills" row or table.

### Checklist for Adding a New Agent

1. Implement the agent (location TBD — likely a new file under `src/datasets_mcp/`).
2. **`PROJECT_INVENTORY.md`** — Add under a new "Agents" section.
3. **`README.md`** — Move from "Planned Agents" to a new "Agents" section. Include an example prompt.
4. **`EXAMPLES.md`** — Add workflow examples showing what the agent does end-to-end.
5. **`ROADMAP.md`** — Move from "Planned Agents" to "Implemented Agents".

### Checklist for Removing / Renaming a Tool

1. Remove or rename in `server.py` (both `list_tools` and `call_tool`).
2. Update the tool set (`BIOPYTHON_TOOLS`, `REST_API_TOOLS`, etc.).
3. Update **all 5 documentation files** — search for the old tool name in each.

## Style Conventions

- Tool names use `snake_case` with a source prefix: `ensembl_`, `uniprot_`, `datasets_`.
- Example prompts in docs are blockquoted (`>`) natural language, followed by a JSON parameter block.
- Parameter tables use: Parameter | Type | Required | Description.
- Section comments in `server.py` use: `# ── Section N: Category ──────────`
- REST tools return curated JSON (extracted fields via `.get()`) — never raw API dumps.
