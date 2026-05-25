# Video prezentacija (5–10 min) – predlog

## Slide 1 (0:00–0:40): Uvod
- Problem: nanoplastici u krvi mogu menjati odgovor imunih ćelija.
- Pitanje: da li veličina čestica (40 vs 200 nm) menja tip i jačinu odgovora?

## Slide 2 (0:40–1:20): Dataset i dizajn
- 4 uzorka istog donora: 40 nm, 200 nm, mix, control.
- scRNA-seq omogućava analizu odgovora po pojedinačnoj ćeliji.

## Slide 3 (1:20–2:30): QC i obrada
- Objasniti filtere i zašto su izabrani.
- Pokazati da su loše ćelije uklonjene pre integracije.

## Slide 4 (2:30–3:30): Integracija i klasteri
- UMAP po uslovu i po klasterima.
- Kratak komentar o batch-correction (Harmony).

## Slide 5 (3:30–4:30): Anotacija cell tipova
- Marker geni + Azimuth PBMC referenca.
- Pomenuti poređenje sa CoDi anotacijama.

## Slide 6 (4:30–5:40): Kompozicija ćelija
- Barplot proporcija cell tipova po uslovima.
- Istaknuti koji tipovi najviše variraju.

## Slide 7 (5:40–7:00): DE i pathways
- Za najveće cell tipove: exposure vs control.
- Fokus na GO/KEGG/Reactome termine sa biološkim smislom.

## Slide 8 (7:00–8:20): Size-specific efekti
- Unique 40, unique 200, shared, mix-only.
- Biološka interpretacija: šta je specifično za veličinu, šta je generalni stress odgovor.

## Slide 9 (8:20–9:20): Dodatne analize
- IFN skor, cell-cycle skor, antigen presentation skor, pseudobulk, agreement metriке.

## Slide 10 (9:20–10:00): Zaključak
- Ključne poruke.
- Ograničenja: jedan donor, potreba za validacijom na više uzoraka.
- Sledeći koraci.
