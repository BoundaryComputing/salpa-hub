# Hello Encrypt

Encrypts a plaintext message using a Caesar cipher and writes the ciphertext to a file.

Part of the [hello-world-pipeline](../README.md) package.

## What it Does

1. Takes a user-provided message (default: "Hello from Salpa!")
2. Applies a Caesar cipher with a configurable shift (default: 3)
3. Writes the encrypted message to `{case_name}_encrypted.txt` with structured markers

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| Case Name | String | `secret` | Identifier for output files; passed to downstream nodes |
| Message | String | `Hello from Salpa!` | The plaintext message to encrypt |
| Cipher Shift | Integer | `3` | Number of positions to shift each letter (1-25) |
| Output Directory | Folder | — | Where to write the encrypted file |
| Include Timestamp | Boolean | `true` | Add timestamp to file header |

## Output

**File**: `{case_name}_encrypted.txt`

```
=== Salpa Secret Message (Encrypted) ===

Case: demo
Cipher: Caesar (shift=3)
Time: 2026-03-26T10:00:00

--- CIPHERTEXT ---
Khoor iurp Vdosd!
--- END ---
```

**Predecessor data** passed to downstream nodes: `{ "case_name": "demo" }`

## Caesar Cipher

Each alphabetic character is shifted forward by `shift` positions in the alphabet, wrapping around from Z to A. Non-alphabetic characters are unchanged.

Example with shift=3: `A→D`, `B→E`, `H→K`, `e→h`, `l→o`
