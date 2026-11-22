#!/usr/bin/env python3
"""Test script for the Datasets MCP Server."""

import asyncio
import json
import os
import sys
import tempfile

import anyio
from mcp.client.session import ClientSession

from datasets_mcp.server import (
    app,
    call_tool as _call_tool,
    check_datasets_installed,
    list_tools as _list_tools,
    run_datasets_command,
)


async def test_basic_functionality() -> int:
    """Run server functionality and protocol tests; return exit code (0=pass, 1=fail)."""
    failures: list[str] = []

    def fail(step: str, msg: str) -> None:
        print(f"   Result: ✗ {msg}")
        failures.append(f"{step}: {msg}")

    print("=" * 60)
    print("Testing Datasets MCP Server")
    print("=" * 60)

    # Test 1: Check if datasets CLI is installed
    print("\n1. Checking if datasets CLI is installed...")
    is_installed = check_datasets_installed()
    if is_installed:
        print("   Result: ✓ Installed")
    else:
        fail("datasets CLI", "datasets CLI not found in PATH")

    if is_installed:
        # Test 2: Get version
        print("\n2. Getting datasets CLI version...")
        result = await run_datasets_command(["version"])
        if result["success"]:
            print(f"   Result: ✓ {result['stdout'].strip()}")
        else:
            fail("datasets version", result["stderr"])

        # Test 3: Test a simple query (summary genome for E. coli)
        print("\n3. Testing genome summary query (E. coli, limit 1)...")
        result = await run_datasets_command([
            "summary", "genome", "taxon", "escherichia coli", "--limit", "1"
        ])
        if result["success"]:
            try:
                data = json.loads(result["stdout"])
                assemblies = data.get("reports", [])
                if assemblies:
                    assembly = assemblies[0]
                    print("   Result: ✓ Found assembly")
                    print(f"      Accession: {assembly.get('accession', 'N/A')}")
                    print(f"      Organism: {assembly.get('organism', {}).get('organism_name', 'N/A')}")
                else:
                    fail("datasets_summary_genome", "Query succeeded but returned no reports")
            except json.JSONDecodeError as e:
                fail("datasets_summary_genome", f"Invalid JSON output: {e}")
        else:
            fail("datasets_summary_genome", result["stderr"])

        # Test 4: Test gene summary query (BRCA1)
        print("\n4. Testing gene summary query (BRCA1 in human)...")
        result = await run_datasets_command([
            "summary", "gene", "symbol", "BRCA1", "--taxon", "human"
        ])
        if result["success"]:
            try:
                data = json.loads(result["stdout"])
                reports = data.get("reports", [])
                if reports:
                    gene = reports[0].get("gene", {})
                    print("   Result: ✓ Found gene")
                    print(f"      Gene ID: {gene.get('gene_id', 'N/A')}")
                    print(f"      Symbol: {gene.get('symbol', 'N/A')}")
                    print(f"      Description: {gene.get('description', 'N/A')[:60]}...")
                else:
                    fail("datasets_summary_gene", "Query succeeded but returned no reports")
            except json.JSONDecodeError as e:
                fail("datasets_summary_gene", f"Invalid JSON output: {e}")
        else:
            fail("datasets_summary_gene", result["stderr"])

        # Test 5: batch_gene_summary
        print("\n5. Testing batch_gene_summary (BRCA1,TP53 in human)...")
        try:
            result = await _call_tool("batch_gene_summary", {
                "symbols": "BRCA1,TP53",
                "taxon": "human",
            })
            output = json.loads(result[0].text)
            assert output["total_results"] >= 2, f"Expected >=2 results, got {output['total_results']}"
            print(f"   Result: ✓ Got {output['total_results']} gene results")
            for r in output["results"].get("reports", []):
                gene = r.get("gene", {})
                print(f"      {gene.get('symbol', 'N/A')}: {gene.get('description', 'N/A')[:50]}")
        except Exception as e:
            fail("batch_gene_summary", str(e))

    # Test 6: sequence_align (nucleotide)
    print("\n6. Testing sequence_align (nucleotide, global)...")
    try:
        result = await _call_tool("sequence_align", {
            "sequence1": "ATCGATCGATCG",
            "sequence2": "ATCGTTCGATCG",
        })
        output = json.loads(result[0].text)
        print(f"   Result: ✓ Alignment score={output['score']}, identity={output['identity_percent']}%")
    except Exception as e:
        fail("sequence_align (nucleotide)", str(e))

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
    except Exception as e:
        fail("sequence_align (protein)", str(e))

    # Test 8: sequence_stats (nucleotide & protein)
    print("\n8. Testing sequence_stats (nucleotide & protein)...")
    try:
        res_nuc = await _call_tool("sequence_stats", {"sequence": "ATCGATCGATCGATCGATCG"})
        out_nuc = json.loads(res_nuc[0].text)
        assert out_nuc["length"] == 20 and out_nuc["gc_content"] == 50.0

        res_prot = await _call_tool("sequence_stats", {
            "sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
            "sequence_type": "protein",
        })
        out_prot = json.loads(res_prot[0].text)
        assert out_prot["length"] == 51
        print(f"   Result: ✓ Nuc GC={out_nuc['gc_content']}%, Prot MW={out_prot.get('molecular_weight', 'N/A')}")
    except Exception as e:
        fail("sequence_stats", str(e))

    # Test 9: parse_fasta
    print("\n9. Testing parse_fasta...")
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1 Human hemoglobin alpha\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH\n")
            f.write(">seq2 Cow hemoglobin alpha\nMVLSAADKGNVKAAWGKVGGHAAEYGAEALERMFLSFPTTKTYFPHFDLSH\n")
            f.write(">seq3 Chicken myoglobin\nMGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLK\n")
            tmp_path = f.name

        try:
            res = await _call_tool("parse_fasta", {"file_path": tmp_path, "include_sequences": False})
            out = json.loads(res[0].text)
            assert out["returned"] == 3

            res2 = await _call_tool("parse_fasta", {"file_path": tmp_path, "filter_pattern": "hemoglobin"})
            out2 = json.loads(res2[0].text)
            assert out2["returned"] == 2
            print(f"   Result: ✓ Parsed {out['returned']} sequences, filter matched {out2['returned']}")
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
    except Exception as e:
        fail("parse_fasta", str(e))

    # Test 10: sequence_translate
    print("\n10. Testing sequence_translate...")
    try:
        coding_dna = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCTGA"
        res = await _call_tool("sequence_translate", {"sequence": coding_dna, "frames": "all"})
        out = json.loads(res[0].text)
        translations = out.get("translations", [])
        assert len(translations) == 6, f"Expected 6 frames, got {len(translations)}"
        assert any(t.get("orfs_10aa_plus") for t in translations), "Expected at least one ORF >= 10 aa"
        print(f"   Result: ✓ Translated {len(translations)} frames")
    except Exception as e:
        fail("sequence_translate", str(e))

    # Test 11: primer_design
    print("\n11. Testing primer_design...")
    try:
        template_500bp = (
            "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA"
            "ATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGAC"
            "CACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAG"
            "TGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGT"
            "CAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAG"
            "TATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGAT"
            "GAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGT"
            "GAACCCGAAAATCCTTCCTTGGAAACCAGTCTCAGTGTCCAACTCTCTAACCTTGGAACT"
            "GTGAGAACTCTGAGGACAAA"
        )
        res = await _call_tool("primer_design", {"sequence": template_500bp})
        out = json.loads(res[0].text)
        assert out.get("primers_found", 0) > 0, f"No primers returned: {out}"
        print(f"   Result: ✓ Designed {out['primers_found']} primer pairs")
    except Exception as e:
        fail("primer_design", str(e))

    # Test 12: hpo_search (all 4 categories)
    print("\n12. Testing hpo_search (search, term, genes, diseases)...")
    try:
        r_search = json.loads((await _call_tool("hpo_search", {"query": "seizure", "category": "search", "limit": 5}))[0].text)
        assert r_search["total_terms"] > 0, f"No HPO terms for 'seizure': {r_search}"

        r_term = json.loads((await _call_tool("hpo_search", {"query": "HP:0001250", "category": "term"}))[0].text)
        assert r_term["id"] == "HP:0001250" and r_term["name"], f"Unexpected term payload: {r_term}"

        r_genes = json.loads((await _call_tool("hpo_search", {"query": "HP:0001250", "category": "genes", "limit": 5}))[0].text)
        assert len(r_genes["genes"]) > 0, f"No HPO genes returned: {r_genes}"

        r_diseases = json.loads((await _call_tool("hpo_search", {"query": "HP:0001250", "category": "diseases", "limit": 5}))[0].text)
        assert len(r_diseases["diseases"]) > 0, f"No HPO diseases returned: {r_diseases}"

        print(
            f"   Result: ✓ search={r_search['total_terms']} terms, "
            f"term='{r_term['name']}', genes={r_genes['total_genes']}, diseases={r_diseases['total_diseases']}"
        )
    except Exception as e:
        fail("hpo_search", str(e))

    # Test 13: pharmgkb_search (ClinPGx)
    print("\n13. Testing pharmgkb_search (gene, drug case-insensitive, variant)...")
    try:
        r_gene = json.loads((await _call_tool("pharmgkb_search", {"gene": "CYP2C9", "limit": 5}))[0].text)
        assert r_gene["total_found"] > 0, f"No ClinPGx results for CYP2C9: {r_gene}"

        r_drug = json.loads((await _call_tool("pharmgkb_search", {"drug": "Warfarin", "limit": 5}))[0].text)
        assert r_drug["total_found"] > 0, f"No ClinPGx results for Warfarin: {r_drug}"

        r_var = json.loads((await _call_tool("pharmgkb_search", {"variant": "rs1799853", "limit": 5}))[0].text)
        assert r_var["total_found"] > 0, f"No ClinPGx results for rs1799853: {r_var}"

        print(
            f"   Result: ✓ CYP2C9={r_gene['total_found']}, "
            f"Warfarin={r_drug['total_found']}, rs1799853={r_var['total_found']}"
        )
    except Exception as e:
        fail("pharmgkb_search", str(e))

    # Test 14: disease_ontology_search (name, synonym, DOID)
    print("\n14. Testing disease_ontology_search (name, synonym, DOID)...")
    try:
        r_name = json.loads((await _call_tool("disease_ontology_search", {"query": "breast cancer", "limit": 5}))[0].text)
        assert r_name["total_results"] > 0, f"No DO results for 'breast cancer': {r_name}"

        r_syn = json.loads((await _call_tool("disease_ontology_search", {"query": "breast tumor", "limit": 5}))[0].text)
        assert r_syn["total_results"] > 0, f"No DO synonym results for 'breast tumor': {r_syn}"

        r_doid = json.loads((await _call_tool("disease_ontology_search", {"query": "DOID:1612"}))[0].text)
        assert r_doid["total_results"] == 1 and r_doid["results"][0]["doid"] == "DOID:1612", f"Bad DOID lookup: {r_doid}"

        print(
            f"   Result: ✓ name={r_name['total_results']}, "
            f"synonym={r_syn['total_results']}, DOID:1612='{r_doid['results'][0]['name']}'"
        )
    except Exception as e:
        fail("disease_ontology_search", str(e))

    # Test 15: Verify tool count directly
    print("\n15. Verifying direct list_tools() count...")
    try:
        tools = await _list_tools()
        assert len(tools) == 45, f"Expected 45 tools, got {len(tools)}"
        print(f"   Result: ✓ Server provides {len(tools)} tools")
    except Exception as e:
        fail("list_tools", str(e))

    # Test 16: End-to-end MCP ClientSession handshake
    print("\n16. Testing MCP ClientSession handshake (initialize, list_tools, call_tool)...")
    try:
        c2s_send, c2s_recv = anyio.create_memory_object_stream(10)
        s2c_send, s2c_recv = anyio.create_memory_object_stream(10)

        async def _run_server() -> None:
            await app.run(c2s_recv, s2c_send, app.create_initialization_options())

        async with anyio.create_task_group() as tg:
            tg.start_soon(_run_server)
            async with ClientSession(s2c_recv, c2s_send) as session:
                await session.initialize()
                mcp_tools = await session.list_tools()
                assert len(mcp_tools.tools) == 45, f"Handshake expected 45 tools, got {len(mcp_tools.tools)}"
                res = await session.call_tool("sequence_stats", {"sequence": "ATCGATCGATCG"})
                payload = json.loads(res.content[0].text)
                assert payload["length"] == 12 and payload["gc_content"] == 50.0
            tg.cancel_scope.cancel()
        print("   Result: ✓ MCP handshake initialized, listed 45 tools, and executed sequence_stats")
    except Exception as e:
        fail("MCP ClientSession handshake", str(e))

    print("\n" + "=" * 60)
    if failures:
        print(f"Test Summary: FAILED ({len(failures)} failure(s))")
        for f_msg in failures:
            print(f"  - {f_msg}")
        print("=" * 60)
        return 1

    print("Test Summary:")
    print("  All 16 tests passed!")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(test_basic_functionality()))
