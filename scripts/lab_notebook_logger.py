#!/usr/bin/env python3
"""
Lab Notebook Logger — Claude Code PostToolUse hook script.

Reads tool call JSON from stdin and appends a timestamped JSONL entry
to .lab-notebook/<session-id>.jsonl in the working directory.

Setup: Add this to your Claude Code hooks configuration:
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
"""

import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# Maximum characters to store from a tool response
MAX_RESPONSE_LENGTH = 2000

# Directory for lab notebook logs
NOTEBOOK_DIR = ".lab-notebook"


def get_session_id() -> str:
    """Get or create a session ID based on the current date and a counter."""
    notebook_dir = Path(NOTEBOOK_DIR)
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    # Find existing sessions for today
    existing = sorted(notebook_dir.glob(f"{today}_*.jsonl"))
    if existing:
        # Use the most recent session file
        return existing[-1].stem
    else:
        return f"{today}_001"


def ensure_notebook_dir() -> Path:
    """Create the .lab-notebook directory if it doesn't exist."""
    notebook_dir = Path(NOTEBOOK_DIR)
    notebook_dir.mkdir(exist_ok=True)
    return notebook_dir


def truncate_response(response: str) -> str:
    """Truncate long responses to keep logs manageable."""
    if len(response) <= MAX_RESPONSE_LENGTH:
        return response
    return response[:MAX_RESPONSE_LENGTH] + f"... [truncated, {len(response)} chars total]"


def main():
    try:
        raw = sys.stdin.read()
        if not raw.strip():
            return

        data = json.loads(raw)
    except (json.JSONDecodeError, Exception):
        return

    notebook_dir = ensure_notebook_dir()
    session_id = get_session_id()
    log_file = notebook_dir / f"{session_id}.jsonl"

    # Extract tool call information from the hook payload
    tool_name = data.get("tool_name", data.get("name", "unknown"))
    tool_input = data.get("tool_input", data.get("arguments", data.get("params", {})))
    tool_response = data.get("tool_response", data.get("response", data.get("result", "")))

    # Truncate response if it's a string
    if isinstance(tool_response, str):
        tool_response = truncate_response(tool_response)
    elif isinstance(tool_response, dict):
        response_str = json.dumps(tool_response)
        if len(response_str) > MAX_RESPONSE_LENGTH:
            tool_response = truncate_response(response_str)

    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "epoch": time.time(),
        "type": "tool_call",
        "tool_name": tool_name,
        "tool_input": tool_input,
        "tool_response": tool_response,
        "session_id": session_id,
    }

    with open(log_file, "a") as f:
        f.write(json.dumps(entry) + "\n")


if __name__ == "__main__":
    main()
