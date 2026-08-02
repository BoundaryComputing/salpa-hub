"""
hello-encrypt — Encrypts a message using a configurable Caesar cipher.

Part of the "Secret Message" hello-world workflow template.
Demonstrates: StringParameter, IntegerParameter, BooleanParameter,
file output, and passing case_name to downstream nodes.

No external dependencies — uses Python stdlib only.
Runs IN_PROCESS (no pixi.toml needed).
"""

import os
from datetime import datetime

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FolderParameter,
    IntegerParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log


def caesar_encrypt(text, shift):
    """Encrypt text using a Caesar cipher with the given shift."""
    result = []
    for char in text:
        if char.isalpha():
            base = ord("A") if char.isupper() else ord("a")
            result.append(chr((ord(char) - base + shift) % 26 + base))
        else:
            result.append(char)
    return "".join(result)


class HelloEncrypt(Node):
    """
    Encrypts a plaintext message using a Caesar cipher.

    Input: Optional predecessor data (uses case_name from upstream if available)
    Output: A text file containing the encrypted message
    """

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            docstring="Name for this case (used in output file names). "
            "If left empty, uses predecessor data or defaults to 'secret'.",
        ),
        "message": StringParameter(
            "Message",
            default="Hello from Salpa!",
            docstring="The plaintext message to encrypt.",
        ),
        "shift": IntegerParameter(
            "Cipher Shift",
            default=3,
            docstring="Number of positions to shift each letter (1-25). "
            "Classic Caesar cipher uses 3.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            default="rel:./",
            docstring="Directory where the encrypted file will be written.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting encryption...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()

            # Read parameters
            case_name = flow_vars["case_name"].get_value()
            message = flow_vars["message"].get_value() or "Hello from Salpa!"
            shift = flow_vars["shift"].get_value()
            if shift is None:
                shift = 3
            shift = max(1, min(25, int(shift)))
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())

            # Use predecessor case_name if ours is empty
            if not case_name and predecessor_data:
                input_data = predecessor_data[0] if predecessor_data else None
                if input_data and isinstance(input_data, dict):
                    case_name = input_data.get("case_name", "secret")
            if not case_name:
                case_name = "secret"

            result.metadata.update({
                "case_name": case_name,
                "execution_time": datetime.now().isoformat(),
            })

            stream_log(f"Encrypting with shift={shift}...", node_id=self.node_id, progress=30)

            # Encrypt the message
            encrypted = caesar_encrypt(message, shift)

            stream_log(f"'{message}' -> '{encrypted}'", node_id=self.node_id, progress=60)

            # Create output directory and write plain ciphertext file
            os.makedirs(output_dir, exist_ok=True)

            output_file = os.path.join(output_dir, f"{case_name}_encrypted.txt")
            with open(output_file, "w") as f:
                f.write(encrypted + "\n")

            stream_log(f"Wrote: {case_name}_encrypted.txt", node_id=self.node_id, progress=90)

            # Build result
            output_path = self.format_output_path(output_file)
            result.data = {
                "case_name": case_name,
                "output_file": output_path,
            }

            result.files["output"] = {
                "encrypted": output_path,
            }

            result.success = True
            result.message = f"Message encrypted for {case_name}"

            stream_log("Encryption complete", node_id=self.node_id, progress=100)
            return result.to_json()

        except Exception as e:
            stream_log(f"Error: {str(e)}", node_id=self.node_id, progress=100)
            raise NodeException("hello-encrypt", str(e))
