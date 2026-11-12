<!-- Author: Jakub Zapletal -->
<!-- Date: 2026-11-12 -->
<!-- Krátký návod, jak pracovat se skripty balíku flood_defenses. -->

# Krátký návod: práce se skripty

Metoda, klíče konfigurace a výstupy jsou popsané v `README.md`. Schémata postupu (mermaid) a vysvětlující obrázky jsou v jeho úvodu, v §5.1 a v §8.

## 1. Příprava

1. **Prostředí pixi** podle `pixi.toml`. Kontrola: `pixi run test-kernels`. Ve výstupu
   musí být jen řádky `PASS`.
2. **Konfigurace:** vyplňuje se jen krátký soubor `config/defenses.toml`:
   - **`[general]`:** `country` (kód země), `work_dir` (rychlý disk pro mezivýsledky),
     `output_dir`, `n_workers` (počet pracovníků × paměť jednotky se musí vejít do RAM);
   - **`[countries.<kód>]`:** cesty k DTM, hrázím a řekám (soubory, které připravíte),
     název sloupce s plochou povodí v km² a volitelně vodicí linie řek (`river_guides`,
     např. vodní toky), podle kterých se řeky přichytí do koryta;
   - **`[countries.<kód>.rp]`:** tabulky škod.

   Ostatní hodnoty jsou ve `config/advanced.toml`. Klíč zapsaný v `defenses.toml` má
   přednost.
   - **Cesty ve Windows:** s lomítky (`"D:/data/dtm.tif"`).
   - **`placeholders`:** klíče s nepotvrzenými hodnotami; skripty je hlásí jako varování.
3. **Kontrola před během:** `pixi run check`. Vypíše OK nebo PROBLEM u každého skriptu a
   vstupu (soubory, sloupce, CRS, složky).

4. **Report po běhu:** `pixi run report` zapíše do repozitáře `reports/<země>/report.md`
   (souhrn vstupů, výstupů a logů, max. 60 KB) a `overview.png`. Složku commitněte a
   pushněte, ať je výsledek vidět i jinde.

## 2. Spouštění

**GUI místo shellu:** `pixi run gui` otevře okno aplikace Defenses:
- **Pipeline:** kroky v pořadí; u každého stav, zapsané výstupy, proměnné běhu (`UNITS`,
  `FORCE`…, výchozí hodnoty ze skriptu) a tlačítko Run. **Run core pipeline** zařadí
  základní kroky za sebou.
- **Konzole** dole: fronta běhů, živý výstup, průběh jednotek s odhadem času, Stop
  (ukončí krok i jeho pracovní procesy).
- **Results:** chráněná území A a B, úniky podle příčin, hráze na koruně, rychlost a odhad
  času celé země, přehled z posledního reportu.
- **Configuration:** formulář pro `defenses.toml` (cesty s kontrolou existence a výběrem
  souboru, sloupec plochy povodí ze sloupců řek); záložka Raw TOML pro celý soubor.
  Komentáře v souboru zůstanou.
- **Logs:** logy skriptů s filtrem varování a chyb a hledáním.

Běh z GUI je stejný jako ze shellu (stejný Python, proměnné `DEFENSES_*`). Zavření:
tlačítko Quit nebo Ctrl+C v terminálu.

**Shell:**

- **Odkud:** vždy z kořene repozitáře: `pixi run <úloha>`.
- **Proměnné běhu** (PowerShell):

  ```powershell
  $env:DEFENSES_COUNTRY = "CZ"          # země
  $env:DEFENSES_UNITS = "012_034,012_035"   # jen vybrané jednotky
  $env:DEFENSES_FORCE = "1"             # přepočítat hotové
  Remove-Item Env:\DEFENSES_FORCE       # zrušení
  ```

  `DEFENSES_FORCE` po použití vždy zrušte.
- **Logy:** `output_dir/<země>/logs/`.
- **Chybná konfigurace:** skript vypíše všechny chyby a nic nepočítá.

## 3. Pořadí

| Krok | `pixi run …` | Kontrola |
|---|---|---|
| 1 | `prepare-dtm` | `dtm_inspection.csv`, `dtm_benchmark.csv` (s DSM i `dsm_*`) |
| 2 | `build-units` | `corridor.gpkg`, `units_core.gpkg` |
| 3 | `protected-areas`: nejdřív dvě jednotky (`DEFENSES_UNITS`), pak vše | `rivers_channel.gpkg` (řeky v korytě), `levees_crest.gpkg` (hráze na koruně), `river_groups.gpkg`, `protected_areas.gpkg` (soubor A), `protected_areas_by_group.gpkg` (soubor B), `leaks.gpkg` (úniky), `barriers.gpkg`, `performance.json` |
| 4 | `download-ghsl`, `rp-index` | `protected_areas_rp.gpkg` (`rp`, `rp_low`, `rp_high`) |
| volitelně | `download-worldcover`, `hand-yz`, `lulc.enabled = true`, znovu `rp-index` | `exposure_lulc.csv` |
| volitelně | `download-gfm`, `gfm-evidence` | `gfm_evidence.gpkg` |
| Stage F | `download-discharge`, `discharge-prepare`, `rp-bounds` | `download_log.csv`, `discharge_cells_<zdroj>.gpkg`, `event_rp.csv`, `obs_rp_low`, `obs_rp_high`, `plots/discharge/base/index.html` |
| volitelně | `validate` | `validation.json` |
| Stage D | `download-clms`, `context-zones` | `corine_codes.tif`, `urban_atlas.gpkg`, `zones.gpkg` |
| Stage E1 | `snap-rivers` | `rivers_snapped.gpkg`, `river_points.gpkg`, `plots/implied/base/index.html` |
| Stage E2 | `implied-segments` | `subsegments.gpkg` |
| Stage E3–E4 | `implied-protection` | `implied/base/implied_subsegments.gpkg`, `implied_protected_areas.gpkg` |
| Stage E5 | `implied-evaluate` | `sensitivity.csv`, `plots/implied/sensitivity/index.html` |

**Soubory A a B** (`protected-areas`, README §8.3):
- **A** (`protected_areas.gpkg`): všechny řeky jsou zdrojem vody najednou.
- **B** (`protected_areas_by_group.gpkg`): chráněná území zvlášť pro každou skupinu toků. Skupina
  jsou celé úseky MERIT spojené po proudu, dokud plocha povodí na soutoku nevzroste víc než
  2× a skupina nepřesáhne 30 km; úseky pod 100 km² tvoří společnou skupinu 0. Hodnoty jsou
  jen v kódu (`gis.Config.DEFAULTS`), ne v toml.
- **Podmínka:** `river_upstream_area_field` musí být **název sloupce** s plochou povodí v km².
  Bez něj se B nepočítá (varování v logu).
- **Čas:** každá skupina v jednotce přidá zhruba tolik, kolik trvá HAND a obě záplavy souboru A
  (na jednotce 30 M buněk asi 27–31 s). Počet skupin je v řádku každé jednotky v logu.
- **Další kroky se souborem B:** `rp-index`, `hand-yz`, `gfm-evidence`, `rp-bounds` a
  `validate` mají proměnnou `LAYER` (`"A"` výchozí, `"B"`). Pro B: nejdřív
  `$env:DEFENSES_LAYER = "B"`, pak tyto kroky ve stejném pořadí; výstupy mají příponu
  `_b` (např. `protected_areas_by_group_rp.gpkg`, `gfm_evidence_b.csv`). V GUI volba
  **File** na kartách těchto kroků.

**Hráze a úniky** (`protected-areas`, README §5.8–5.10):
- **`levees_crest.gpkg`:** linie hrází posunuté na korunu z DTM (nejvýš 15 m). Sloupec
  `no_ridge` označí linie, které v DTM jako hráz nejsou vidět (chyba mapy, násep, který
  není hráz, nebo hráz pod rozlišením DTM).
- **Mezery:** mezery mezi linkami do 200 m se spojí jen tam, kde DTM drží korunu.
  V `barriers.gpkg` mají `bridge_type = "dtm"`.
- **`leaks.gpkg`:** u hrází, které nic nechrání, místo a příčina:
  - `river_crossing`: otevřený průchod toku hrází;
  - `low_crest`: voda přeteče nízkou korunou;
  - `line_end`: voda obteče volný konec linie;
  - `water_source_behind`: tok uvnitř, který hrází neprochází;
  - `terrain_above_crest`: terén za hrází je výš než koruna (hráz tam není potřeba,
    nebo je koruna podhodnocená).

  `levee_length_m` je délka hráze, která tím místem uniká, `lost_margin_median_m` ochrana,
  kterou únik bere. Shrnutí podle příčin a pět největších úniků je v logu. Začněte
  největšími úniky.
- **Po změně linií hrází** spusťte `protected-areas` s FORCE (stejně jako u řek).

**Náspy silnic a železnic** (`transport-embankments`, README §5.11, zatím jen analýza):
- **Vstup:** do `[countries.<kód>]` přidejte `transport_lines = ["D:/…/roads.gpkg"]`, tedy
  linie silnic a železnic, které připravíte. Ukazují jen, kde hledat; o všem rozhoduje DTM.
- **Kdy:** spouští se po `protected-areas` (potřebuje HAND jednotek pro nivu).
- **Výstupy:**
  - `transport_embankments.gpkg`: třídy násep, jednostranný, zářez, v úrovni terénu;
  - `openings.gpkg`: mosty, podjezdy, pravděpodobné propustky s jistotou;
  - grafy v `plots/transport/<běh>/index.html`.
- **Parametry:** jsou v `[transport]` v `config/advanced.toml`; ověřte je na pilotních
  jednotkách podle řezů a podélných profilů.
- **Bariéry:** náspy do bariér zatím nevstupují. To bude další krok (integrace).

**Stage E** je experimentální a nejdřív se pouští na pilotní jednotky (`[implied] units`). Výsledky se do
`defenses_clip.gpkg` nedostanou, dokud není `include_implied = true`.

**Stage F** (průtoky GloFAS a EFAS z EWDS):
- **Přístup.** Osobní klíč patří do `~/.cdsapirc` nebo do proměnné `CDSAPI_KEY`, nikdy do
  repozitáře. Licence obou datasetů je potřeba přijmout ve formuláři EWDS.
- **Hodnoty požadavků** (`request`, `data_format`, roky, `time_keys`, `area_key`) se
  kopírují z formuláře EWDS („Show API request“). Prázdné položky skript ohlásí a zdroj
  přeskočí.
- **Kontrola časové konvence.** V logu `discharge-prepare` je řádek „valid times … ->
  days …“. Ověřte ho na staženém souboru a na známé povodni.
- **Opakované spuštění** `download-discharge` stáhne jen chybějící nebo neúplné soubory.
- **Bez EFAS.** Když EFAS není dostupný, běží vše jen s GloFAS.

## 4. Po změně spustit znovu

| Změna | Spustit |
|---|---|
| hráze, řeky, vodicí linie, parametry chráněných území | `protected-areas` s FORCE, pak `rp-index` (a `hand-yz`, `gfm-evidence`, `rp-bounds`) |
| tabulky škod, výška podlaží | `rp-index` s FORCE, pak `rp-bounds` |
| třídy, strop, váhy, ruční zásahy | `rp-index`, pak `rp-bounds` |
| hodnoty land cover | `hand-yz`, pak `rp-index`, `rp-bounds` |
| nové satelitní scény | `download-gfm`, `gfm-evidence`, `rp-bounds` |
| parametry `[discharge]`, nové roky | `download-discharge`, `discharge-prepare`, `rp-bounds` |
| řeky (Stage E) | `snap-rivers`, `implied-segments` s FORCE, `implied-protection` |
| parametry `[implied]`, budovy | `implied-protection`, pak `implied-evaluate` |

Po každém `rp-index` znovu `rp-bounds`, pokud se používají pozorování.

## 5. Časté problémy

- **`contains N layers`:** vektor má více vrstev. Uložte vrstvu samostatně (QGIS export,
  `ogr2ogr`).
- **Nedostatek paměti:** snižte `n_workers` nebo `unit_core_size_m`, pak znovu
  `build-units`.
- **Neověřené URL:** data stáhněte ručně.
  - **GHSL:** šablona `file://…`.
  - **GloFAS, EFAS:** soubory z formuláře EWDS uložte do `<složka>/<zdroj>/<rok>[_<měsíc>].<přípona>` a
    nastavte `ewds_api_url = "file:///<složka>"` (ve Windows např. `file:///C:/data/ewds`).
- **Pomalé čtení DTM:** `convert_dtm = true` (varianta A), pak znovu `prepare-dtm`.
- **Kontrola celého řetězce:** `test-pipeline`, `test-stage-b`, `test-stage-c`,
  `test-hand-yz`, `test-stage-d`, `test-stage-e`, `test-stage-f`.
- **Grafy Stage E a F:** `index.html` ve složce `plots/implied/<běh>/`, resp.
  `plots/discharge/<běh>/`. Zobrazení během výpočtu: `$env:DEFENSES_SHOW_PLOTS = "1"`.
- **Pomalé stahování EWDS.** Rychlost určuje fronta služby. Měsíc 6hodinových dat EFAS
  může trvat dlouho; `max_parallel_requests` zvyšujte opatrně.
