# Reed–Solomon Error‑Correcting Code 

A complete Python implementation of Reed–Solomon encoding and decoding (single‑error correction, s=1). This project demonstrates how to:

- Encode any text message into a Reed–Solomon codeword (blocks of 20 characters → integers mod p)
- Detect and correct a single bit‑flip error in the codeword
- Decode and recover the original message

---

## Features

- RS encode/decode over a large prime field (161‑bit prime)
- Three methods for computing the free coefficient (naive, k‑optimized, one‑inversion)
- Automatic error detection & correction via Lagrange interpolation
- Polynomial arithmetic utilities (addition, multiplication, scaling)
- Conversion between text ↔ integers ↔ bytes

---
