# Plan rada

Ovaj dokument definira plan rada sekvenciranja gena iz FASTQ datoteka.

---

## 1. DIO: Priprema i postavljanje
Prikupljanje potrebnih resursa i vanjskih biblioteka.
- [x] **Preuzimanje podataka:** Osigurati testne `.fastq` datoteke.
- [x] **Instalacija biblioteka:**
    - [x] [spoa](https://github.com/rvaser/spoa) (Simulated Partial Order Alignment)
    - [x] [bioparser](https://github.com/rvaser/bioparser) (Fastq/Fasta parser)

---

## 2. DIO: Upravljanje podatcima (C++)
Implementacija memorijskih struktura za efikasno čitanje genoma.
1. **Struktura `Sequence`**:
   - Ova će struktura sadržavati polja iz FASTQ datoteka koja će parse pročitati (Name, Data, Quality).
2. **Batch Loading**:
   - Čitanje u blokovima radi uštede RAM-a (1 GB po iteraciji kao u primjeru).
   - Svi se objekti pokranjuju u vektor `std::vector<std::unique_ptr<Sequence>>`.

---

## 3. DIO: Analiza i usporedba 
Glavni dio logike koji izračunava sličnost između parova sekvenci.

### Klasa za Analizu Sekvenci
Razred će uspoređivati svaku sekvencu sa svakom (izuzev sa samom sobom) kako bi utvrdio najmanju udaljenost.

| Zadatak | Opis |
| :--- | :--- |
| **Pronađi Min Udaljenost** | Za svaku sekvencu $S_i$ pronađi sekvencu $S_j$ s najmanjim Edit Distance-om tako da $i\ne j$. |
| **Mapiranje Rezultata** | Metoda koja vraća mapirani ID na uređeni par (tuple). Tuple se sastoji od ID-a odgovarajućeg para sekvence za koji je pronađena najmanja udaljenost te iznos najmanje udaljenosti. |
| **Struktura Izlaza** | `std::map<ID, std::tuple<Neighbor_ID, Distance>>` |

---

## 4. DIO: Clustering
Završna obrada podataka i grupiranje.
- [ ] Implementirati računanje **centroida** (reprezentativna sekvenca grupe).
- [ ] Implementirati **clustering** algoritam temeljen na prethodno izračunatim udaljenostima.

---
