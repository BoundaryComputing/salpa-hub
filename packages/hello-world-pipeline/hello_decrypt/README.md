# Hello Decrypt

Reads a Caesar-cipher encrypted file and writes the decrypted plaintext to a new file.

Part of the [hello-world-pipeline](../README.md) package.

## What it Does

1. Receives `case_name` from the upstream node (Hello Encrypt)
2. Reads `{case_name}_encrypted.txt` and extracts the ciphertext between markers
3. Reverses the Caesar cipher with the configured shift
4. Writes the plaintext to `{case_name}_decrypted.txt`

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| Case Name | String | *(from predecessor)* | Used to locate the encrypted file |
| Cipher Shift | Integer | `3` | Must match the encryption shift to correctly decode |
| Input Directory | Folder | — | Where to find the encrypted file |
| Output Directory | Folder | — | Where to write the decrypted file |
| Include Timestamp | Boolean | `true` | Add timestamp to file header |

## Output

**File**: `{case_name}_decrypted.txt`

```
=== Salpa Secret Message (Decrypted) ===

Case: demo
Cipher: Caesar (shift=3)
Time: 2026-03-26T10:00:01

--- PLAINTEXT ---
Hello from Salpa!
--- END ---
```

**Predecessor data** passed to downstream nodes: `{ "case_name": "demo" }`

## Error Handling

- If the encrypted file is not found, the node raises a clear error with the expected file path and a hint to check that Hello Encrypt ran first
- If the ciphertext markers (`--- CIPHERTEXT ---` / `--- END ---`) are missing, a descriptive error is raised
