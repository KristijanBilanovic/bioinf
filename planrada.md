# Plan rada
- preuzeti podatke
- preuzeti ```spoa```    
- preuzeti ```bioparser``` 
- implementirati čitanje sekvence gena u C++ strukturu
- spremiti sve pročitane sekvence u vektor
- napaviti klasu koja uspoređuje svaku sekvencu sa svakom, pamti najmanje udaljenosti
  - ova klasa mora imati metodu koja vraća mapirani ID svake pročitane sekvence na tuple (uređeni par) koji se sastoji od ID-ja druge sekvence (s kojom je pronađena minimalna udaljenosti)s te iznosa udaljenosti
- implementirati računanje centroiida i clusterings