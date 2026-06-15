# Skript za odbranu — šta da kažeš za svaku sliku i tabelu

**Projekat:** Single-cell analiza imunog odgovora na nanoplastične čestice  
**Podaci:** 33,240 ćelija, 4 uslova (control, 40 nm, 200 nm, mix)  
**Napomena:** Brojevi su iz run-a `pipeline_run_20260611_212256`. Jedan donor — uvek dodaj „u ovom uzorku“ ili „hipoteza za dalja istraživanja“.

---

# DEO A — SLIKE (figures)

---

## 1. `umap_condition.png` — UMAP po eksperimentalnom uslovu

### Šta je ovo (1 rečenica za publiku)
„UMAP je mapa gde je svaka tačka jedna ćelija, a blizina znači sličan genetski profil. Ovde su ćelije obojene po tome da li su bile na 40 nm, 200 nm, mix-u ili control-u.“

### Šta da pokažeš prstom
- Sve četiri boje su **rasute po istim regionima** — nema jedne boje koja zauzima ceo ugao ili odvojeno ostrvo.
- Postoje „oblaci“ (T ćelije, B ćelije, monociti) — to je **normalna PBMC struktura**, ne uslov.

### Šta konkretno da kažeš
„Vidite da se control, 40 nm, 200 nm i mix **preklapaju** na istim delovima mape. To znači da nanoplastika **nije potpuno promenila identitet** svake ćelije na globalnom nivou — mapa izgleda slično za sve uslove. To **nije loš rezultat**: integracija je uspešna, i biološki je očekivano da in vitro izloženost ne razbije celu populaciju krvnih ćelija. Ali to **ne znači da nema efekta** — efekti su vidljivi u diferencijalnoj ekspresiji gen po gen, po tipu ćelije, što pokazujemo u tabelama.“

### Šta je očekivano za krv (PBMC)
- U perifernoj krvi **dominiraju T ćelije** (~40–50%) — to ćemo videti i na sledećim slikama.
- B ćelije, NK i monociti formiraju **odvojene regije** na UMAP-u — tipičan PBMC obrazac.

### Ako te pitaju
**„Zašto se ne vide odvojeno?“** → Efekat je suptilan na nivou celog transcriptoma; jaki efekti su u specifičnim genima (DE), ne u pomeranju cele ćelije na mapi.

---

## 2. `umap_split_by_condition.png` — UMAP podeljen po uslovima (4 panela)

### Šta je ovo
„Ista mapa, ali u četiri panela: u svakom panelu jedan uslov obojen, ostalo sivo u pozadini.“

### Šta da kažeš
„Ovaj prikaz je **lakši za oko** nego jedna preklopljena slika. U svakom panelu vidite da obojene tačke **pokrivaju iste oblasti** koje su sive u pozadini — npr. T ćelije, B ćeoka, monociti. Nijedan uslov nema ekskluzivan region samo za sebe. To potvrđuje: **struktura populacije ostaje**, nanoplastika ne pravi potpuno novu populaciju ćelija.“

### Detalj po uslovu (možeš kratko)
| Uslov | Šta da primetiš |
|-------|-----------------|
| Control | Referentni raspored — svi tipovi prisutni |
| 40 nm | Isti layout; nema velikog pomeranja |
| 200 nm | Isti layout; proporcije se malo menjaju (videti barplot) |
| Mix | Isti layout; mix ponekad daje drugačije gene, ne drugačiju mapu |

---

## 3. `umap_sample_integration.png` — UMAP po uzorku (batch)

### Šta je ovo
„Svaka boja = jedan tehnički uzorak (Sample 1–4). Pre Combat korekcije uzorci mogu biti odvojeni; posle Combat treba da se **mešaju**.“

### Šta da kažeš
„Combat je uklonio tehničke razlike između četiri uzorka. Ako bismo videli četiri odvojena ostrva, **ne bismo mogli pouzdano** da poredimo uslove — mešali bismo batch sa biologijom. Ovde se uzorci **preklapaju**, što znači da integracija radi i da UMAP po uslovu (prethodna slika) zaista pokazuje biologiju, a ne artefakt sekvenciranja.“

### Očekivano
- Dobar integrisan dataset → **nema** jasne granice po `sample_id`.
- Loš integrisan → svaki sample u svom uglu (to **nemamo**).

---

## 4. `umap_clusters.png` — Leiden klasteri (14 klastera)

### Šta je ovo
„Klasteri su **nekontrolisane grupe** — algoritam grupiše slične ćelije bez znanja o tipu. Dobili smo **14 klastera**.“

### Šta da kažeš
„Pre nego što smo dodelili T/B/NK/monocit, algoritam je sam našao 14 grupa. Broj je **razuman za PBMC** — nije previše sitno (šum), nije premalo (gubitak detalja). Klasteri odgovaraju glavnim kompartmanima krvi. To govori da je **kvalitet podataka dobar** — struktura u podacima postoji i nije samo šum.“

### Očekivano za PBMC
- 10–20 klastera na ~30k ćelija sa resolution 0.5 je uobičajeno.
- Kasnije se klasteri povezuju sa marker genima.

---

## 5. `umap_celltypes_marker.png` — Tipovi ćelija (marker anotacija)

### Šta je ovo
„Svaka ćelija dobila tip na osnovu **kanonskih marker gena** iz literature (CD3 za T, MS4A1 za B, LYZ za monocite…).“

### Brojevi koje **moraš** da znaš (ceo dataset, svi uslovi)

| Tip | Broj ćelija | % |
|-----|------------:|--:|
| CD4 T | 15,568 | **46.8%** |
| NK | 4,064 | 12.2% |
| B | 3,986 | 12.0% |
| CD8 citotoksične T | 3,449 | 10.4% |
| Monocit CD14 | 2,313 | 7.0% |
| DC | 2,078 | 6.3% |
| Monocit CD16 | 1,488 | 4.5% |
| Trombociti | 294 | 0.9% |

### Šta da kažeš
„**CD4 T ćelije su najbrojnije** — oko 47% svih ćelija. To je **potpuno očekivano za PBMC** iz periferne krvi: T limfociti su dominantna populacija. B ćelije su u **odvojenom ostrvu** — klasičan znak dobre anotacije. Monociti su u drugom delu mape. **NK i CD8 T se delimično preklapaju** — to je normalno jer delimo gene poput NKG7 u marker panelu; nisu potpuno odvojive populacije na UMAP-u.“

### Očekivano za krv vs naše
- **Očekivano u PBMC:** CD4 T > CD8 T; B ~5–15%; NK ~5–15%; monociti ~10–20% (zavisi od preparata).
- **Kod nas:** profil odgovara tipičnom PBMC — nema npr. 80% monocita (što bi bilo sumnjivo).

### Napomena za DC
„Imamo 6.3% DC po markerima — **više nego Azimuth** (46 ćelija na L1). Marker panel za DC (`FCER1A`, `CST3`) može **preceniti** DC ili uhvatiti ćelije na granici monocit–DC. To spominjem kao **ograničenje**, ne kao grešku u podacima.“

---

## 6. `umap_codi_celltypes.png` — CoDi spoljna anotacija

### Šta je ovo
„CoDi su **nezavisne etikete** iz Zenodo CSV fajlova — druga metoda, ne naši markeri. Mapirano **99.6%** ćelija.“

### Šta da kažeš
„Raspored na mapi je **sličan marker anotaciji** — T regioni, B ostrvo, monociti. To je **validacija**: dve metode se slažu u glavnim linijama. Ne moraju biti identične ćelija po ćelija — CoDi i markeri imaju drugačiju granularnost.“

### Broj
- **99.6%** ćelija ima CoDi labelu — odlično mapiranje.

---

## 7. `umap_module_scores.png` — Skorovi (S, G2M, IFN)

### Šta je ovo
„Umesto kategorije, boja = **broj (skor)** koliko je ćelija u tom programu: proliferacija (S, G2M) ili interferon (IFN).“

### Šta da kažeš po panelu

**S_score i G2M_score (cell cycle):**  
„Tamnije/viša vrednost = više gena faze deljenja. Vidite da su ti signali **lokalizovani** u određenim regionima — ne ravnomerno po celoj mapi. To je normalno: samo mali deo PBMC u krvi aktivno deli. **Globalno po uslovima** razlike su male (videti tabelu cell_cycle).“

**IFN_score:**  
„Interferon potpis — innate imun odgovor. Prostor na mapi pokazuje gde IFN program postoji; **nije** da ceo UMAP „gore“ pod nanoplastikom. U tabeli vidimo da je control **malo viši** od izloženih — ne pričam priču o masivnoj IFN aktivaciji od PSNP.“

### Povezivanje sa tabelama
- Cell-cycle CSV: S i G2M blizu control-a.
- IFN CSV: control 0.035, 200 nm 0.010, 40 nm 0.005, mix −0.002.

---

## 8. `marker_dotplot.png` — Dot plot marker gena

### Šta je ovo
„Veličina tačke = koliko ćelija tipa **ekspresira** gen; boja = **koliko** (intenzitet).“

### Šta da kažeš
„Ovo je **dokaz da anotacija ima smisla**. Na primer: `MS4A1`, `CD79A` su visoki u B ćelijama; `LYZ`, `S100A8` u CD14 monocitima; `NKG7` u NK i delimično CD8 T. **T_cell markeri** (`CD3D`, `CD3E`) u T populacijama. Bez ovoga komisija može pitati ‚odakle znaš da su to B ćelije?‘ — ovde je odgovor.“

### Očekivano za PBMC markere
| Gen / tip | Očekivanje |
|-----------|------------|
| MS4A1, CD79A | B ćelije |
| CD3D, CD3E | T ćelije |
| LYZ, S100A8 | Monociti |
| NKG7, GNLY | NK / citotoksični T |
| PPBP, PF4 | Trombociti |

Ako dot plot to poštuje → **anotacija je validna**.

---

## 9. `composition_barplot.png` — Sastav populacije po uslovu

### Šta je ovo
„Za svaki uslov — **udeo** svakog tipa ćelije (100% = jedan stubac).“

### Glavna priča — **porasti i padovi** (control → izložen)

| Tip | Control | 40 nm | 200 nm | Mix | Trend |
|-----|--------:|------:|-------:|----:|-------|
| **CD4 T** | **49.2%** | 48.0% | **44.0%** ↓ | 48.4% | Pad kod **200 nm** (~5 pp) |
| **Monocit CD14** | **5.7%** | 6.6% | **9.7%** ↑↑ | **3.1%** ↓ | **Najveći porast kod 200 nm** (+4 pp) |
| NK | 12.7% | 12.0% | 11.3% | **13.9%** | Blage promene |
| B | 11.3% | 12.1% | 11.9% | 12.8% | Stabilno |
| CD8 T | 10.6% | 9.8% | 10.6% | 10.5% | Stabilno |
| Monocit CD16 | 3.7% | 4.0% | **5.2%** | 4.5% | Blagi porast 200 nm |
| DC | 5.5% | 7.2% | 5.9% | 6.4% | Umerene oscilacije |
| Trombociti | 1.3% | 0.2% | 1.4% | 0.3% | **Mali broj ćelija** — ne preciziraj jako |

### Šta da kažeš (detaljno)
„**Globalno, sastav PBMC ostaje sličan control-u** — i dalje dominiraju T ćelije, što je očekivano za krv. **Ključna promena:** kod **200 nm** PSNP udeo **CD14 monocita raste sa 5.7% na 9.7%** — skoro duplo u relativnom smislu. To sugerise da **veće čestice (200 nm)** možda jače angažuju **mieloidnu granu** (innate, fagocitoza). Kod **40 nm** monociti su blago viši (6.6%), ne tako dramatično. **Mix** ima **najmanje** CD14 monocita (3.1%) — zanimljivo, mix ne kopira 200 nm; **mešavina veličina može imati drugačiji efekat** nego pojedinačne veličine. CD4 T padaju na 200 nm (49 → 44%) — može biti kompenzacija povećanjem monocita u istom stubcu.“

### Očekivano za krv
- T uvek najveći udeo — **kod nas svuda 44–49%** ✓  
- Monociti obično **10–20%** u „klasičnom“ PBMC; kod nas **3–10%** zavisno od uslova — preparat i donor utiču; **relativni** skok na 200 nm je važniji od apsolutne vrednosti.

### Oprez
„Jedan donor — kažem **u ovom donoru** / **hipoteza**, ne univerzalni zakon.“

---

# DEO B — TABELE (tables)

---

## 1. `cell_composition_by_condition.csv`

**Ista priča kao barplot**, ali sa tačnim brojevima ćelija.

### Primeri brojeva koje možeš citirati
- Control CD4 T: **3,129** ćelija (49.2%)
- 200 nm CD4 T: **5,470** ćelija (44.0%) — više ćelija u apsolutnom broju jer je sample veći, ali **udeo** pada
- Control CD14 mono: **364** (5.7%) → 200 nm: **1,205** (9.7%)

### Kako da objasniš tabelu
„Svaki red = jedan uslov × jedan tip. Kolona `fraction` = udeo u tom uslovu. Koristimo ovo da **kvantifikujemo** ono što vidimo na barplot-u. Glavni zaključak: **proporcije su uglavnom stabilne**, izuzev **porasta CD14 monocita na 200 nm**.“

---

## 2. `differential_expression_all.csv`

### Šta je ovo
„**66,000 redova** — za svaku ćeliju tipa i poređenje (40/200/mix vs control) testirani genovi: log fold change, p-vrednost, adjusted p.“

### Brojevi
- **22 poređenja** (8 tipova × do 3 uslova; trombociti samo 200 nm vs control)
- **36,814 značajnih** gena (padj < 0.05, |logFC| > 0.25)

### Šta da kažeš
„Dok UMAP izgleda slično, **na nivou gena efekat je masivan**: skoro 37 hiljada značajnih promena. To znači da nanoplastika **menja ekspresiju** — ali **ne isto u svim tipovima** i **ne isti genovi** za 40 vs 200 nm. Wilcoxon test unutar svakog tipa — poredimo 40 nm T sa control T, ne T sa B. Metoda: **Wilcoxon**, praga **padj < 0.05**, **|logFC| > 0.25**.“

### Očekivano / interpretacija
- Jak DE + stabilan UMAP = **specifične promene u programima**, ne potpuna promena identiteta ćelije.
- **Jedan donor** — p-vrednosti su informativne unutar dataset-a, ali treba više donora za generalizaciju.

### Primer kako da zvučiš stručno
„DE je **stratifikovana po tipu ćelije** — to je ispravan pristup za PBMC jer T i monocit ne smeju biti mešani u jednom testu.“

---

## 3. `pathway_enrichment_all.csv`

### Šta je ovo
„**89,015 redova** — za značajne DE gene, Enrichr testira da li su previše zastupljeni u GO, KEGG, Reactome putanjama.“

### Šta da kažeš
„Putanje koje se **najčešće** pojavljuju govore o **imunologiji i inflamaciji**. Konkretni primeri iz naših rezultata:“

| Tip / uslov | Primer putanje |
|-------------|----------------|
| NK, 200 nm vs control | **Inflammatory Response**, **Cytokine-mediated signaling** |
| B, 200 nm vs control | **Inflammatory Response** |
| CD8 T, 200 nm | **Cytokine-mediated signaling** |
| Monocit CD14, **mix** vs control | **Inflammatory Response**, **Cytokine-mediated signaling** |
| Monocit CD14, 200 nm | Regulation of cell migration, translation (različit profil od mix-a) |

„To **podržava priču** da nanoplastika aktivira **imune gene programe** — citokini, inflamacija — posebno kod **200 nm i mix-a**, i u **NK, B, CD8, monocitima**. To je **očekivano** ako čestice interaguju sa imunim ćelijama u krvi.“

### Očekivano za imun izazov
- Inflamatorni odgovor, hemotaksija, citokini — **plauzibilno** za strani materijal / NP.
- Nije dokaz in vivo patologije — **in vitro** PBMC.

---

## 4. `size_specific_effects_summary.csv`

### Šta je ovo
„Broji DE gene po **klasi efekta**: samo 40 nm, samo 200 nm, deljeni 40+200, deljeni sva tri, **samo mix** (emergentni).“

### Najvažniji redovi — **citiraj ove brojeve**

| Tip | unique 40nm | unique 200nm | shared 40+200 | **shared all 3** | mix only |
|-----|------------:|-------------:|--------------:|-----------------:|---------:|
| **Monocit CD14** | 203 | 81 | 267 | **865** | 294 |
| **NK** | 205 | 287 | 398 | **748** | 253 |
| Monocit CD16 | 221 | 269 | 487 | 697 | 183 |
| CD8 T | 180 | 328 | 310 | 527 | **382** |
| DC | 163 | 296 | 177 | 504 | **402** |
| B | 198 | 335 | 339 | 624 | 270 |
| CD4 T | 325 | 260 | 470 | 452 | 261 |
| Trombociti | 0 | **1419** | 0 | 0 | 0 |

### Šta da kažeš (detaljno)
„**Veličina čestice važi:** na primer B ćelije imaju **335 gena samo za 200 nm** i **198 samo za 40 nm** — to nisu isti genovi. **Monocit CD14** ima **865 gena** značajnih **u sva tri** izloženja vs control — to tretiram kao **jezgro odgovora** klasičnih monocita na PSNP. **Mix** dodaje **294 gena** kod CD14 koje **nema** ni 40 nm ni 200 nm solo — **emergentni efekat mešavine**. DC i CD8 T imaju **400+ mix-only gena** — mešavina nije samo prosek veličina. **Trombociti:** skoro sve je **unique 200 nm** (1419 gena) — imamo **samo jedno** DE poređenje za platelets, pa tu **ne generalizujem** o veličini.“

### Poruka za odbranu
1. **Size-specific** — da, 40 ≠ 200 nm.  
2. **Shared core** — monociti/NK nose najveći zajednički set.  
3. **Mix emergent** — mix nije trivijalan.

---

## 5. `cell_cycle_scores_by_condition.csv`

### Šta je ovo
„Prosečan **S** i **G2M** skor po uslovu (proliferacija / deljenje).“

### Brojevi

| Uslov | S_score | G2M_score |
|-------|--------:|----------:|
| control | 0.0029 | 0.0120 |
| 40 nm | 0.0049 | 0.0163 |
| 200 nm | 0.0005 | **−0.0177** |
| mix | 0.0044 | 0.0122 |

### Šta da kažeš
„Skorovi su **blizu nule** i **mali** — većina PBMC u krvi **nije u aktivnom deljenju**, što je **normalno** (T i B u mirnom stanju). **Nema dramatične** promene ciklusa kod nanoplastike: 40 nm i mix blago slični control-u; 200 nm ima **malo niži G2M** (−0.018 vs 0.012). To možda sugerise blagu **supresiju proliferacije** na 200 nm, ali razlike su **male** — ne gradim glavni zaključak na ovome. Važno: u ovom run-u skorovi **rade** (nisu NaN) jer su računati pre HVG filtriranja.“

### Očekivano za krv
- Periferna krv: **nizak** proliferativni indeks vs tkivo — **kod nas niski skorovi** ✓

---

## 6. `ifn_scores_by_condition.csv`

### Brojevi

| Uslov | IFN_score |
|-------|----------:|
| **control** | **0.0351** (najviši) |
| 200 nm | 0.0096 |
| 40 nm | 0.0045 |
| mix | −0.0018 |

### Šta da kažeš
„IFN potpis meri innate / antiviral program (`ISG15`, `IFIT1`, `MX1`…). **Control je nešto viši** od izloženih — suprotno od jednostavne priče ‚nanoplastika pali interferon‘. Razlike su **slabe**. Kažem **iskreno**: IFN modul **nije glavni dokaz** efekta PSNP; jači su DE i pathway rezultati. Ne prikrivam negativan nalaz.“

### Očekivano
- Jak IFN od NP → izloženi > control; **kod nas nije tako** → **slab / neodlučan** modul.

---

## 7. `antigen_presentation_scores.csv`

### Šta je ovo
„HLA / antigen prezentacija — **MHC put**; po uslovu i tipu ćelije.“

### Očekivano gde su skorovi visoki
- **B ćelije** i **profesionalni APC** (monociti, DC) — MHC program.

### Konkretni brojevi (B ćelije — najviši)

| Uslov | B_cell skor |
|-------|------------:|
| control | **1.44** |
| mix | 1.32 |
| 200 nm | 1.17 |
| 40 nm | 0.95 |

Monocit CD14 control **0.26** → mix **0.12** (pad); CD16 control **0.43** → 40 nm **0.29**.

### Šta da kažeš
„**B ćelije imaju najviše** antigen-prezentacione skorove — **očekivano** (MHC II, CD74). Control B je na **1.44**, 40 nm pada na **0.95** — moguća **blaga supresija** MHC programa na 40 nm kod B. Monociti i DC imaju **niže apsolutne vrednosti** ali su biološki relevantni za **nanoplastiku i imunitet** — gledamo **relativne** promene po tipu, ne samo globalni prosek. Mix monocit CD14 **0.12** vs control **0.26** — pad HLA signala u mix-u za klasične monocite.“

---

## 8. `annotation_agreement_metrics.csv`

### Broj
- **64.7%** slažnost CoDi vs marker tip

### Šta da kažeš
„Dve **nezavisne** metode se slažu u **roughly dve trećine** ćelija. **Ne očekujem 100%** — CoDi kaže ‚CD4+ T cell‘, mi ‚CD4_T‘; NK i citotoksični T dele markere. **64.7% je u redu** za validaciju: glavne linije (T, B, mono) se poklapaju, fine granice ne moraju.“

---

## 9. `pseudobulk_counts_condition_celltype.csv`

### Šta je ovo
„**32 grupe** — sabrani UMI brojevi po (uslov × tip ćelije). Kao ‚bulk‘ sažetak za spoljnu validaciju ili model.“

### Šta da kažeš
„Tehnički output za **downstream** analizu — npr. DESeq2 na pseudobulk-u sa više donora u budućnosti. **Ne prikazujem na slajdu** osim ako pitaju; spomenem da postoji za replikaciju studije.“

---

## 10. `azimuth_annotations.csv` (Azimuth)

### Brojevi
- **33,240** ćelija  
- Mean score **0.908**, median **0.967**  
- Samo **250** ćelija (<0.75%) sa score < 0.5  

### L1 raspodela (top)
| Tip | Ćelije |
|-----|-------:|
| CD8 T | 13,108 |
| CD4 T | 11,568 |
| B | 2,652 |
| NK | 2,382 |
| Mono | 1,789 |
| DC (L1) | **46** |

### Šta da kažeš
„Azimuth mapira na **referentni PBMC atlas** — mean confidence **0.91** je **odlično**. Potvrđuje T dominaciju i B/NK/mono. **DC samo 46** na L1 vs **2,078** po markerima — ista priča kao na UMAP-u: marker DC panel **precenjuje** u odnosu na referencu. Koristim Azimuth za **finije T podtipove** i **pouzdanost**, ne kao jedini izvor DC broja.“

---

# DEO C — REDOSLED ZA 10-MINUTNU ODBRANU

1. **Dizajn** — 1 donor, 4 uslova, 33k ćelija, Zenodo  
2. **umap_condition + split** — struktura ista, nema ostrva  
3. **umap_sample_integration** — Combat OK  
4. **umap_celltypes + marker_dotplot** — PBMC profil očekivan, CD4 ~47%  
5. **composition_barplot** — **monociti ↑ 200 nm (5.7→9.7%)**  
6. **DE** — 36,814 hitova  
7. **pathway + size_specific** — inflamacija, veličina važi, mix emergent  
8. **cell_cycle / IFN** — iskreno: cycle stabilan, IFN slab  
9. **Azimuth + CoDi** — validacija 0.91 / 64.7%  
10. **Limitacije** — jedan donor, in vitro  

---

# DEO D — ENGLISH PHRASES (ako odbrana na engleskom)

| Srpski koncept | English one-liner |
|----------------|-------------------|
| UMAP preklapanje | „All conditions occupy the same UMAP territories — no global transcriptome collapse.“ |
| Monocit porast | „CD14 monocytes increase from 5.7% to 9.7% under 200 nm PSNP in this donor.“ |
| DE vs UMAP | „Strong DE despite similar UMAP — effects are gene- and cell-type-specific.“ |
| Size | „865 shared DE genes in CD14 monocytes across all exposures — a core response module.“ |
| IFN | „IFN module scores are weak and slightly higher in control — not our main line of evidence.“ |
| Limitation | „Single-donor exploratory analysis; findings require replication.“ |

---

*Kraj skripta. Prilagodi ton profesoru — brojeve citiraj sa slajda ili tabele pored tebe.*
