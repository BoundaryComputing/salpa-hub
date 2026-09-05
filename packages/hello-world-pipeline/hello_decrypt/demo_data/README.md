# demo_data

`demo/demo_encrypted.txt` — *Khoor iurp Vdosd!*, written by the `hello_encrypt` node (message
*Hello from Salpa!*, shift 3, case name `demo`), so it is exactly what this node reads from
its input directory upstream. `DEMO_CONFIG` in `node.py` points the node at `demo_data/demo`
with the same case name and shift. Decrypting it gives *Hello from Salpa!* back.
