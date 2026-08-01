# Hello World Pipeline

A 3-node workflow that encrypts, decrypts, and reveals a secret message using a Caesar cipher. This is the simplest multi-node pipeline in Salpa — designed to teach you how workflows work without any external dependencies.

## What You'll Learn

- **Connecting nodes** — wire output ports to input ports on the canvas
- **Data flow** — how `case_name` propagates from node to node via predecessor data
- **File-based communication** — each node reads/writes files in a shared working directory
- **Node configuration** — each node has its own independently configurable parameters
- **Sequential execution** — nodes run in topological order, waiting for predecessors

## The Pipeline

```
[Hello Encrypt] ──→ [Hello Decrypt] ──→ [Hello Reveal]
```

| Step | Node | What it does |
|------|------|-------------|
| 1 | **Hello Encrypt** | Takes your message, encrypts it with a Caesar cipher, writes `{case}_encrypted.txt` |
| 2 | **Hello Decrypt** | Reads the encrypted file, reverses the cipher, writes `{case}_decrypted.txt` |
| 3 | **Hello Reveal** | Reads both files, generates a visual comparison report in `{case}_reveal.txt` |

## Quick Start

### Option A: Load from Template (recommended)

1. Open the **Load Workflow** dialog in the Salpa canvas
2. Switch to the **Templates** tab
3. Find **"Hello World Pipeline"** and click **Load**
4. All three nodes are pre-wired and configured — just click **Execute All**

The bundled workflow template (`workflows/hello-world-pipeline.json`) sets up the full pipeline with sensible defaults.

### Option B: Build Manually

1. **Create a new workflow** in the Salpa canvas
2. **Add the three nodes** from the node palette (search "hello")
3. **Connect them** in order: Encrypt → Decrypt → Reveal
4. **Configure** each node:
   - Set all `Output Directory` / `Input Directory` fields to `rel:./` (the workflow's working folder)
   - On the Encrypt node, set your `Message` (e.g., "Hello from Salpa!") and `Cipher Shift` (default: 3)
   - On the Decrypt node, set the same `Cipher Shift` value
   - Set a `Case Name` on the first node (e.g., "demo") — it will flow to downstream nodes automatically
5. **Execute** the workflow — click Execute on each node in order, or use Execute All

## Output Files

After execution, your working directory will contain:

```
demo_encrypted.txt    ← Ciphertext with file markers
demo_decrypted.txt    ← Recovered plaintext
demo_reveal.txt       ← Visual comparison report
```

### Example Report (detailed style)

```
==================================================
  SECRET MESSAGE — DETAILED ANALYSIS
==================================================

  Case:    demo
  Cipher:  Caesar (shift=3)

--------------------------------------------------
  ENCRYPTED:
    Khoor iurp Vdosd!

  DECRYPTED:
    Hello from Salpa!

--------------------------------------------------
  CHARACTER MAPPING:

    Encrypted:  Khoor iurp Vdosd!
                ||||||||||||||||||||
    Decrypted:  Hello from Salpa!

--------------------------------------------------
  STATISTICS:
    Total characters:   20
    Alphabetic:         16
    Shifted:            16
    Unchanged:          4

  Status: Successfully decoded
==================================================
```

## Node Parameters

### Hello Encrypt

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| Case Name | String | `secret` | Identifier for output files; passed to downstream nodes |
| Message | String | `Hello from Salpa!` | The plaintext message to encrypt |
| Cipher Shift | Integer | `3` | Caesar cipher shift (1-25) |
| Output Directory | Folder | — | Where to write the encrypted file |
| Include Timestamp | Boolean | `true` | Add timestamp to file header |

### Hello Decrypt

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| Case Name | String | *(from predecessor)* | Used to locate the encrypted file |
| Cipher Shift | Integer | `3` | Must match the encryption shift |
| Input Directory | Folder | — | Where to find the encrypted file |
| Output Directory | Folder | — | Where to write the decrypted file |
| Include Timestamp | Boolean | `true` | Add timestamp to file header |

### Hello Reveal

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| Case Name | String | *(from predecessor)* | Used to locate both files |
| Input Directory | Folder | — | Where to find encrypted + decrypted files |
| Output Directory | Folder | — | Where to write the report |
| Report Style | Select | `detailed` | `box` (compact) or `detailed` (with character mapping) |
| Include Timestamp | Boolean | `true` | Add timestamp to report |

## How Data Flows

Only `case_name` is passed between nodes via Salpa's predecessor data mechanism. All other parameters are configured independently on each node.

```
Encrypt                    Decrypt                    Reveal
┌──────────────────┐      ┌──────────────────┐      ┌──────────────────┐
│ case_name: "demo"│─────→│ case_name: "demo"│─────→│ case_name: "demo"│
│ message: "Hello" │      │ shift: 3         │      │ report: detailed │
│ shift: 3         │      │ input: rel:./    │      │ input: rel:./    │
│ output: rel:./   │      │ output: rel:./   │      │ output: rel:./   │
└──────────────────┘      └──────────────────┘      └──────────────────┘
        │                         │                         │
        ▼                         ▼                         ▼
  demo_encrypted.txt       demo_decrypted.txt        demo_reveal.txt
```

## Technical Details

- **No external dependencies** — pure Python stdlib
- **Execution strategy**: IN_PROCESS (no pixi.toml, no environment isolation needed)
- **Cross-platform**: Works on Linux, macOS, and Windows
- **Caesar cipher**: A simple substitution cipher that shifts each letter by a fixed number of positions in the alphabet. Non-alphabetic characters (spaces, punctuation, numbers) are preserved unchanged.

## Package Info

- **Package**: `hello-world-pipeline`
- **Version**: 1.0.0
- **License**: MIT
- **Author**: Boundary Computing
