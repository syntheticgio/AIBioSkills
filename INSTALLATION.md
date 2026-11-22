# Installation & Setup Guide

## Prerequisites

1. **Python 3.10+** — with `pip` or [`uv`](https://docs.astral.sh/uv/)
2. **NCBI Datasets CLI** *(optional for BioPython and REST API tools; required for the NCBI Datasets CLI tools)* — install following the [NCBI Datasets CLI installation docs](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)

## Installation Steps

### 1. Install the MCP Server

Clone the repository and install in editable mode using `pip` or `uv`:

#### Option A: Using `pip` and `venv`

```bash
git clone https://github.com/syntheticgio/AIBioSkills.git
cd AIBioSkills
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -e .
```

#### Option B: Using `uv`

```bash
git clone https://github.com/syntheticgio/AIBioSkills.git
cd AIBioSkills
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .
```

### 2. Optional Environment Variables

- **`SEMANTIC_SCHOLAR_API_KEY`** *(optional)*: Provide a Semantic Scholar API key to obtain higher rate limits when using `semantic_scholar_search`. Without a key, public rate limits apply.

  ```bash
  export SEMANTIC_SCHOLAR_API_KEY="your-api-key"
  ```

### 3. Verify Installation

Test that everything works (with your virtual environment activated):

```bash
# Verify imports
python -c "from datasets_mcp.server import app; print('OK')"

# Run the test script
python test_server.py

# Check that the entry point is on PATH
which datasets-mcp-server
```

Expected output from `test_server.py`:
```
✓ datasets CLI installed (<version>)
✓ Server provides 45 tools
✓ Genome queries working (tested with E. coli)
✓ Gene queries working (tested with BRCA1)
```

### 4. Configure Claude Desktop

#### Linux Configuration

Edit `~/.config/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "datasets": {
      "command": "datasets-mcp-server"
    }
  }
}
```

**Note:** If `datasets-mcp-server` is not in your global `PATH` (for example, when installed inside a project virtual environment), use the full path to the executable inside `.venv`:

```json
{
  "mcpServers": {
    "datasets": {
      "command": "/path/to/AIBioSkills/.venv/bin/datasets-mcp-server"
    }
  }
}
```

To pass the optional `SEMANTIC_SCHOLAR_API_KEY` for higher Semantic Scholar rate limits, add an `"env"` block:

```json
{
  "mcpServers": {
    "datasets": {
      "command": "/path/to/AIBioSkills/.venv/bin/datasets-mcp-server",
      "env": {
        "SEMANTIC_SCHOLAR_API_KEY": "your-api-key"
      }
    }
  }
}
```

#### macOS Configuration

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` using the same JSON configuration as above.

### 5. Restart Claude Desktop

After editing the configuration:
1. Completely quit Claude Desktop (not just close the window)
2. Restart Claude Desktop
3. The `datasets` MCP server should now be available

## Verification

Once configured, you can test in Claude Desktop by asking:

- "Can you check the datasets tool version?"
- "Get information about the BRCA1 gene in humans"
- "Show me genome assemblies for E. coli"

## Troubleshooting

### Server not appearing in Claude Desktop

1. Check the config file location and syntax (must be valid JSON)
2. Verify the command path is correct
3. Check Claude Desktop logs for errors
4. Try the absolute path to `datasets-mcp-server` (e.g., `/path/to/AIBioSkills/.venv/bin/datasets-mcp-server`)

### "datasets CLI not found" error

The NCBI Datasets CLI tools check for `datasets` in `PATH`. Verify:

```bash
which datasets
```

If not found, install it following the [NCBI Datasets CLI installation docs](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) and make sure its directory is on your `PATH`. (BioPython and REST API tools do not require the `datasets` CLI.)

### Import errors

If you get module import errors, reinstall in your active virtual environment from the repository root:

```bash
pip install -e .
# Or with uv:
uv pip install -e .
```

### Network/timeout errors

NCBI Datasets and external bioinformatics APIs require internet access. Check:
- Internet connection is active
- No firewall is blocking NCBI or external API endpoints
- Increase timeout in `src/datasets_mcp/server.py` if needed

## Uninstallation

To remove the MCP server:

```bash
pip uninstall datasets-mcp-server
# Or with uv:
uv pip uninstall datasets-mcp-server
```

Then remove the configuration from `claude_desktop_config.json`.

## Development

To modify the server:

1. Edit files in `src/datasets_mcp/`
2. Changes take effect immediately (editable install with `-e`)
3. Restart Claude Desktop to pick up changes
4. Run `python test_server.py` to verify changes

## Additional Tools

### MCP Inspector

Test the server interactively:

```bash
npx @modelcontextprotocol/inspector datasets-mcp-server
```

This provides a web interface to test tool calls.

### Direct CLI Testing

Test `datasets` commands directly:

```bash
# Check version
datasets version

# Test a query
datasets summary genome taxon human --limit 1

# Test gene query
datasets summary gene symbol BRCA1 --taxon human
```

## Support

For issues:
- Server code issues: Check `README.md` in the repository root
- NCBI Datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
- MCP protocol: https://modelcontextprotocol.io/
