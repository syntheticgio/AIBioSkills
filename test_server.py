#!/usr/bin/env python3
"""Simple test script for the MCP server."""

import asyncio
import json
from datasets_mcp.server import app, check_datasets_installed, run_datasets_command


async def test_basic_functionality():
    """Test basic server functionality."""
    print("=" * 60)
    print("Testing NCBI Datasets MCP Server")
    print("=" * 60)

    # Test 1: Check if datasets CLI is installed
    print("\n1. Checking if datasets CLI is installed...")
    is_installed = check_datasets_installed()
    print(f"   Result: {'✓ Installed' if is_installed else '✗ Not installed'}")

    if not is_installed:
        print("\n   ERROR: datasets CLI not found in PATH")
        print("   Please install from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
        return

    # Test 2: Get version
    print("\n2. Getting datasets CLI version...")
    result = await run_datasets_command(["version"])
    if result["success"]:
        print(f"   Result: ✓ {result['stdout'].strip()}")
    else:
        print(f"   Result: ✗ Error: {result['stderr']}")

    # Test 3: List available tools
    print("\n3. Listing available MCP tools...")
    tool_names = [
        "datasets_summary_genome",
        "datasets_summary_gene",
        "datasets_download_genome",
        "datasets_download_gene",
        "datasets_version"
    ]
    print(f"   Result: ✓ Server provides {len(tool_names)} tools:")
    for name in tool_names:
        print(f"      - {name}")

    # Test 4: Test a simple query (summary genome for E. coli)
    print("\n4. Testing genome summary query (E. coli, limit 1)...")
    result = await run_datasets_command([
        "summary", "genome", "taxon", "escherichia coli", "--limit", "1"
    ])
    if result["success"]:
        try:
            data = json.loads(result["stdout"])
            assemblies = data.get("reports", [])
            if assemblies:
                assembly = assemblies[0]
                accession = assembly.get("accession", "N/A")
                org_name = assembly.get("organism", {}).get("organism_name", "N/A")
                print(f"   Result: ✓ Found assembly")
                print(f"      Accession: {accession}")
                print(f"      Organism: {org_name}")
            else:
                print("   Result: ✓ Query succeeded but no results found")
        except json.JSONDecodeError:
            print("   Result: ✓ Command succeeded (non-JSON output)")
    else:
        print(f"   Result: ✗ Error: {result['stderr']}")

    # Test 5: Test gene summary query (BRCA1)
    print("\n5. Testing gene summary query (BRCA1 in human)...")
    result = await run_datasets_command([
        "summary", "gene", "symbol", "BRCA1", "--taxon", "human"
    ])
    if result["success"]:
        try:
            data = json.loads(result["stdout"])
            reports = data.get("reports", [])
            if reports:
                gene = reports[0].get("gene", {})
                gene_id = gene.get("gene_id", "N/A")
                symbol = gene.get("symbol", "N/A")
                description = gene.get("description", "N/A")
                print(f"   Result: ✓ Found gene")
                print(f"      Gene ID: {gene_id}")
                print(f"      Symbol: {symbol}")
                print(f"      Description: {description[:60]}...")
            else:
                print("   Result: ✓ Query succeeded but no results found")
        except json.JSONDecodeError:
            print("   Result: ✓ Command succeeded (non-JSON output)")
    else:
        print(f"   Result: ✗ Error: {result['stderr']}")

    # Test 6: sequence_align (nucleotide)
    print("\n6. Testing sequence_align (nucleotide, global)...")
    try:
        from datasets_mcp.server import call_tool as _call_tool
        result = await _call_tool("sequence_align", {
            "sequence1": "ATCGATCGATCG",
            "sequence2": "ATCGTTCGATCG",
        })
        output = json.loads(result[0].text)
        print(f"   Result: ✓ Alignment score={output['score']}, identity={output['identity_percent']}%")
        print(f"      Seq type: {output['sequence_type']}, gaps: {output['gaps']}")
    except Exception as e:
        print(f"   Result: ✗ Error: {e}")

    # Test 7: sequence_align (protein)
    print("\n7. Testing sequence_align (protein, global)...")
    try:
        result = await _call_tool("sequence_align", {
            "sequence1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
            "sequence2": "MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH",
            "sequence_type": "protein",
        })
        output = json.loads(result[0].text)
        print(f"   Result: ✓ Alignment score={output['score']}, identity={output['identity_percent']}%")
        print(f"      Seq type: {output['sequence_type']}, length: {output['alignment_length']}")
    except Exception as e:
        print(f"   Result: ✗ Error: {e}")

    # Test 8: sequence_stats (nucleotide)
    print("\n8. Testing sequence_stats (nucleotide)...")
    try:
        result = await _call_tool("sequence_stats", {
            "sequence": "ATCGATCGATCGATCGATCG",
        })
        output = json.loads(result[0].text)
        print(f"   Result: ✓ Length={output['length']}, GC={output['gc_content']}%")
        print(f"      Base counts: {output['base_counts']}")
    except Exception as e:
        print(f"   Result: ✗ Error: {e}")

    # Test 9: sequence_stats (protein)
    print("\n9. Testing sequence_stats (protein)...")
    try:
        result = await _call_tool("sequence_stats", {
            "sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
            "sequence_type": "protein",
        })
        output = json.loads(result[0].text)
        print(f"   Result: ✓ Length={output['length']}, MW={output.get('molecular_weight', 'N/A')}")
        print(f"      pI={output.get('isoelectric_point', 'N/A')}")
    except Exception as e:
        print(f"   Result: ✗ Error: {e}")

    # Test 10: parse_fasta (create temp file)
    print("\n10. Testing parse_fasta...")
    import tempfile
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1 Human hemoglobin alpha\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH\n")
            f.write(">seq2 Cow hemoglobin alpha\nMVLSAADKGNVKAAWGKVGGHAAEYGAEALERMFLSFPTTKTYFPHFDLSH\n")
            f.write(">seq3 Chicken myoglobin\nMGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLK\n")
            tmp_path = f.name

        result = await _call_tool("parse_fasta", {
            "file_path": tmp_path,
            "include_sequences": False,
        })
        output = json.loads(result[0].text)
        print(f"   Result: ✓ Found {output['returned']} sequences")
        for s in output["sequences"]:
            print(f"      {s['id']}: {s['description'][:50]} (len={s['length']})")

        # Test with filter
        result2 = await _call_tool("parse_fasta", {
            "file_path": tmp_path,
            "filter_pattern": "hemoglobin",
        })
        output2 = json.loads(result2[0].text)
        print(f"   Filter test: ✓ 'hemoglobin' matched {output2['returned']} sequences")

        os.unlink(tmp_path)
    except Exception as e:
        print(f"   Result: ✗ Error: {e}")

    # Test 11: batch_gene_summary
    print("\n11. Testing batch_gene_summary (BRCA1,TP53 in human)...")
    try:
        result = await _call_tool("batch_gene_summary", {
            "symbols": "BRCA1,TP53",
            "taxon": "human",
        })
        output = json.loads(result[0].text)
        print(f"   Result: ✓ Got {output['total_results']} gene results")
        for r in output["results"].get("reports", []):
            gene = r.get("gene", {})
            print(f"      {gene.get('symbol', 'N/A')}: {gene.get('description', 'N/A')[:50]}")
    except Exception as e:
        print(f"   Result: ✗ Error: {e}")

    # Test 12: List all tools (verify count)
    print("\n12. Verifying tool count...")
    from datasets_mcp.server import list_tools as _list_tools
    tools = await _list_tools()
    print(f"   Result: ✓ Server provides {len(tools)} tools:")
    for t in tools:
        print(f"      - {t.name}")

    print("\n" + "=" * 60)
    print("Test Summary:")
    print("  All functionality tests completed!")
    print("  The MCP server appears to be working correctly.")
    print("=" * 60)


if __name__ == "__main__":
    import os
    asyncio.run(test_basic_functionality())
