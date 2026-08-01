# Hello Reveal

Reads both the encrypted and decrypted files and generates a visual comparison report.

Part of the [hello-world-pipeline](../README.md) package.

## What it Does

1. Receives `case_name` from the upstream node (Hello Decrypt)
2. Reads both `{case_name}_encrypted.txt` and `{case_name}_decrypted.txt`
3. Extracts and compares the ciphertext and plaintext
4. Generates a formatted report in one of two styles
5. Writes the report to `{case_name}_reveal.txt`

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| Case Name | String | *(from predecessor)* | Used to locate the encrypted/decrypted files |
| Input Directory | Folder | — | Where to find both input files |
| Output Directory | Folder | — | Where to write the report |
| Report Style | Select | `detailed` | `box` for compact summary, `detailed` for character mapping |
| Include Timestamp | Boolean | `true` | Add timestamp to report |

## Report Styles

### `detailed` (default)

Includes character-by-character mapping and statistics:

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

### `box`

Compact framed summary:

```
+==========================================+
|      SECRET MESSAGE REVEALED             |
+------------------------------------------+
|  Case:       demo                        |
|  Cipher:     Caesar (shift=3)            |
|  Encrypted:  Khoor iurp Vdosd!           |
|  Decrypted:  Hello from Salpa!           |
|  Status:     Match verified              |
+==========================================+
```
