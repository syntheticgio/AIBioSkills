---
name: lab-notebook
description: Research session lab notebook — start, annotate, update, report, or check status
user-invocable: true
---

You are a lab notebook assistant. Based on the argument provided (`$ARGUMENTS`), perform one of these actions:

## Subcommands

### `start <title>`
Start a new lab notebook session.
1. Create the `.lab-notebook/` directory if it doesn't exist
2. Determine the session ID: today's date + incrementing counter (e.g., `2026-02-12_002`)
3. Write a new JSONL file `.lab-notebook/<session-id>.jsonl` with an initial entry:
   ```json
   {"timestamp": "...", "type": "session_start", "title": "<title>", "objective": "<title>", "session_id": "<session-id>"}
   ```
4. Confirm the session has started and remind the user to enable the PostToolUse hook if they haven't already:
   ```
   Hook config for .claude/settings.json:
   {
     "hooks": {
       "PostToolUse": [
         {
           "hooks": [
             {
               "type": "command",
               "command": "python3 scripts/lab_notebook_logger.py"
             }
           ]
         }
       ]
     }
   }
   ```

### `annotate <note>`
Add a manual annotation to the current session log.
1. Use the `lab_notebook_annotate` MCP tool with the provided note text
2. Confirm the annotation was added

### `update` or `update <note>`
Read the JSONL log for the current session and produce a narrative summary of recent activity.
1. Read the most recent `.lab-notebook/*.jsonl` file
2. Group tool calls into logical research steps (by time proximity and tool relationships)
3. Write a brief narrative update summarizing what was done, what was found
4. If a note was provided, incorporate it as context
5. Output the update to the user

### `report`
Generate a polished final lab notebook as a markdown file.
1. Read the most recent `.lab-notebook/*.jsonl` file
2. Parse all entries chronologically
3. Group tool calls into logical research steps:
   - Consecutive calls to related tools (e.g., gene lookup → variant search → protein features) form one step
   - Time gaps > 5 minutes suggest a new step
   - Identify **abandoned branches**: sequences of tool calls whose results weren't built upon in subsequent steps — format these with strikethrough
4. Generate a structured markdown notebook with:
   - **Title and Objective** (from `session_start` entry)
   - **Timeline** with numbered steps, each containing:
     - Step title and timestamp
     - Narrative description of what was done and why
     - Tool calls with key results (bullet points)
   - **Abandoned branches** (strikethrough sections with explanation)
   - **Key Findings** (numbered list of discoveries)
   - **Tools Used** (table with tool name and call count)
   - **Session Metadata** (duration, total calls, data sources)
5. Write the report to `lab-notebook-<session-id>.md` in the working directory
6. Display the report to the user

### `status`
Show current session statistics.
1. Read the most recent `.lab-notebook/*.jsonl` file
2. Count total tool calls, unique tools, session duration
3. List the tools used with counts
4. Show the session start time and title
5. Display a concise status summary

## Notes
- The JSONL log is created by the `scripts/lab_notebook_logger.py` hook script
- Each line is a JSON object with: timestamp, type, tool_name, tool_input, tool_response, session_id
- `session_start` entries have: timestamp, type, title, objective, session_id
- `annotation` entries have: timestamp, type, note, session_id
- Always use the most recent session file (sorted by filename) unless the user specifies otherwise
