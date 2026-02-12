# Installation & Setup Guide

## Prerequisites

1. **NCBI Datasets CLI** - Already installed at `/home/syntheticgio/anaconda3/bin/datasets`
   - Version: 18.5.3 ✓

2. **Python 3.10+** - Available through Anaconda ✓

## Installation Steps

### 1. Install the MCP Server

The server is already installed in your conda base environment:

```bash
cd /home/syntheticgio/Programming/datasets-mcp
conda run -n base pip install -e .
```

### 2. Verify Installation

Test that everything works:

```bash
# Run the test script
conda run -n base python test_server.py

# Test the CLI command
datasets-mcp-server --help
```

Expected output from test script:
```
✓ datasets CLI installed (version 18.5.3)
✓ MCP server provides 5 tools
✓ Genome queries working (tested with E. coli)
✓ Gene queries working (tested with BRCA1)
```

### 3. Configure Claude Desktop

#### Linux Configuration

Edit: `~/.config/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "datasets": {
      "command": "datasets-mcp-server"
    }
  }
}
```

**Note:** If `datasets-mcp-server` is not in your PATH, use the full path:

```json
{
  "mcpServers": {
    "datasets": {
      "command": "/home/syntheticgio/anaconda3/bin/datasets-mcp-server"
    }
  }
}
```

Or run via conda:

```json
{
  "mcpServers": {
    "datasets": {
      "command": "conda",
      "args": ["run", "-n", "base", "datasets-mcp-server"]
    }
  }
}
```

#### macOS Configuration

Edit: `~/Library/Application Support/Claude/claude_desktop_config.json`

Use the same JSON configuration as above.

### 4. Restart Claude Desktop

After editing the configuration:
1. Completely quit Claude Desktop (not just close the window)
2. Restart Claude Desktop
3. The datasets MCP server should now be available

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
4. Try the absolute path to datasets-mcp-server

### "datasets CLI not found" error

The server checks for `datasets` in PATH. Verify:

```bash
which datasets
# Should output: /home/syntheticgio/anaconda3/bin/datasets
```

If not found, make sure your conda environment is activated when running the server.

### Import errors

If you get module import errors:

```bash
# Reinstall in the correct environment
cd /home/syntheticgio/Programming/datasets-mcp
conda run -n base pip install -e .
```

### Network/timeout errors

NCBI datasets requires internet access. Check:
- Internet connection is active
- No firewall blocking NCBI APIs
- Try increasing timeout in server.py if needed

## Uninstallation

To remove the MCP server:

```bash
conda run -n base pip uninstall datasets-mcp-server
```

Then remove the configuration from `claude_desktop_config.json`.

## Development

To modify the server:

1. Edit files in `/home/syntheticgio/Programming/datasets-mcp/src/datasets_mcp/`
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

Test datasets commands directly:

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
- Server code issues: Check `/home/syntheticgio/Programming/datasets-mcp/README.md`
- NCBI datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/
- MCP protocol: https://modelcontextprotocol.io/
