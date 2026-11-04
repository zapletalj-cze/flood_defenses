<!-- Author: Jakub Zapletal -->
<!-- Date: 2026-11-04 -->
<!-- Krátký návod, jak pracovat se skripty balíku flood_defenses. -->

# Krátký návod: práce se skripty

Metoda, klíče konfigurace a výstupy jsou popsané v `README.md`.

## 1. Příprava

1. **Prostředí pixi** podle `pixi.toml`. Kontrola: `pixi run test-kernels`. Ve výstupu
   musí být jen řádky `PASS`.
2. **Konfigurace `config/defenses.toml`:**
   - **`[general]`:** `work_dir` (rychlý disk pro mezivýsledky), `output_dir`,
     `n_workers` (počet pracovníků × paměť jednotky se musí vejít do RAM).
   - **Tabulky země:** zkopírujte tabulky `[countries.XX…]`, přepište `XX` na kód země a
     vyplňte cesty.
   - **Cesty ve Windows:** s lomítky (`"D:/data/dtm.tif"`) nebo v apostrofech
     (`'D:\data\dtm.tif'`).
   - **`placeholders`:** klíče s nepotvrzenými hodnotami; skripty je hlásí jako varování.

## 2. Spouštění

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
| 3 | `download-osm`, `osm-hydraulics` | `rivers_sources.gpkg` (`is_source`, `reason`) |
| 4 | `compartments`: nejdřív dvě jednotky (`DEFENSES_UNITS`), pak vše | `compartments.gpkg`, `barriers.gpkg`, `performance.json` |
| 5 | `download-ghsl`, `rp-index` | `compartments_rp.gpkg` (`rp`, `rp_low`, `rp_high`) |
| volitelně | `download-worldcover`, `hand-yz`, `lulc.enabled = true`, znovu `rp-index` | `exposure_lulc.csv` |
| volitelně | `download-gfm`, `gfm-evidence`, `download-glofas`, `rp-bounds` | `gfm_evidence.gpkg`, `obs_rp_low`, `obs_rp_high` |
| volitelně | `validate` | `validation.json` |
| Stage D | `download-clms`, `context-zones` | `corine_codes.tif`, `urban_atlas.gpkg`, `zones.gpkg` |

## 4. Po změně spustit znovu

| Změna | Spustit |
|---|---|
| hráze, řeky, parametry kompartmentů, OSM | `compartments` s FORCE, pak `rp-index` (a `hand-yz`, `gfm-evidence`, `rp-bounds`) |
| tabulky škod, výška podlaží | `rp-index` s FORCE, pak `rp-bounds` |
| třídy, strop, váhy, ruční zásahy | `rp-index`, pak `rp-bounds` |
| hodnoty land cover | `hand-yz`, pak `rp-index`, `rp-bounds` |
| nové satelitní scény | `download-gfm`, `gfm-evidence`, `rp-bounds` |

Po každém `rp-index` znovu `rp-bounds`, pokud se používají pozorování.

## 5. Časté problémy

- **`contains N layers`:** vektor má více vrstev. Uložte vrstvu samostatně (QGIS export,
  `ogr2ogr`).
- **Nedostatek paměti:** snižte `n_workers` nebo `unit_core_size_m`, pak znovu
  `build-units`.
- **Neověřené URL:** data stáhněte ručně.
  - **OSM:** `osm_pbf`.
  - **GHSL:** šablona `file://…`.
  - **GloFAS:** `glofas_files`.
- **Pomalé čtení DTM:** `convert_dtm = true` (varianta A), pak znovu `prepare-dtm`.
- **Kontrola celého řetězce:** `test-pipeline`, `test-stage-b`, `test-stage-c`,
  `test-hand-yz`, `test-stage-d`.
