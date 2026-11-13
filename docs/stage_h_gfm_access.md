<!-- Author: Jakub Zapletal -->
<!-- Date: 2026-11-13 -->
<!-- Stage H, step 1: GFM access and products as verified on 2026-11-13. -->

# Stage H, step 1: GFM access and products

Checked on 2026-11-13 against the live service. Everything below was read from the service
or its Product User Manual (PUM, GFM D5.2, V20251201); what could not be verified is marked.

## Access

| Item | Verified value |
|---|---|
| STAC API | `https://stac.eodc.eu/api/v1` (search by POST `/search`, paged by the `next` link) |
| Collection | `GFM` (one of 58 collections; no separate "scientific" collection is listed) |
| Authentication | none for the STAC search and for the files (`https://data.eodc.eu/...`, read with `/vsicurl/`) |
| Coverage | 2015-01-01 to now; items of 2016, 2020, 2022 and 2024 all carry `processing:version = GFM v3.1.0` (a reprocessed archive) |
| Item | one Sentinel-1 scene on one Equi7 tile of 300 km (15 000 x 15 000 px at 20 m) |
| Item properties | `datetime`, `Equi7Tile`, `flooded_pixels`, `floodable_pixels`, `anomaly_detected`, `flood_members`, `proj:wkt2`, `proj:transform`, `proj:shape` |
| REST API | `https://api.gfm.eodc.eu/v2` (swagger at `/v2/swagger.json`); every data endpoint answers 401 without a bearer token (`/auth/get_bearer_token`, registered users) |
| Portal | `https://portal.gfm.eodc.eu` (registration) |

## Licence

- The STAC collection gives `license: other` and links to the CEMS terms of service
  (`https://global-flood.emergency.copernicus.eu/terms-of-service/`, last modified 2023-06-07).
- Terms: free, full and open access under Regulation (EU) 2021/696; reproduction,
  distribution, communication to the public, adaptation and combination are allowed.
- Required notices:
  - unchanged data: "Generated using Copernicus Emergency Management Service information [Year]";
  - adapted data: "Contains modified Copernicus Emergency Management Service information [Year]".
- No warranty; some data are restricted to registered users (article 53). The files of the
  `GFM` collection were readable without registration.

## Assets of an item

| Asset key | Layer | Values found |
|---|---|---|
| `ensemble_flood_extent` | Observed Flood Extent (final) | 0 no flood, 1 flood, 255 no data |
| `exclusion_mask` | Exclusion Mask (final) | 0, 1 excluded, 255 no data |
| `reference_water_mask` | Reference Water Mask (monthly) | 0, 1, 2 (permanent and seasonal water; which is which is not stated in the PUM), no 255 in the tile read |
| `ensemble_likelihood` | Likelihood (final) | 0–100; flood pixels 50–100, non-flood 0–49 (PUM); 255 no data |
| `advisory_flags` | Advisory Flags | 0, 1 low regional backscatter, 2 rough water, 3 both (PUM); 255 no data |
| `ensemble_water_extent` | Observed Water Extent | flood plus reference water |
| `dlr_*`, `list_*`, `tuw_*` | intermediate results of the three algorithms | not needed |

Checks on the most flooded scene of the pilot event (2024-09-28, tile E051N018T3, a
6000 x 6000 px window):
- flood pixels never coincide with the exclusion mask or the reference water (the product
  sets them to 0), so "wet" must be read as flood only, and reference water is neither wet
  nor a valid dry observation;
- **65.9 % of the valid pixels are excluded** (forest, built-up areas, topography, radar
  shadow); valid observations inside protected areas will be much scarcer than the swath
  suggests;
- likelihood of flood pixels: p5 48, median 94, p95 100.

## CRS

- No EPSG code in the files or in the items: the CRS is the Equi7 grid of Europe, an
  azimuthal equidistant projection given only as WKT
  (`+proj=aeqd +lat_0=53 +lon_0=24 +x_0=5837287.81977 +y_0=2121415.69617 +datum=WGS84`).
- Other continents use other Equi7 zones. The CRS must therefore be read from every file
  (GDAL does this); nothing may assume EPSG:3035 or EPSG:4326.

## File layout and read speed

- Flood extent, exclusion mask, likelihood: COG, 512 px blocks, ZSTD, 5 overviews. A
  6000 x 6000 px window took 0.9–1.8 s.
- Reference water mask: **strip layout (15000 x 1), LZW, no overviews**: the same window
  took 2.9 s, and a small window still reads whole rows. It is one file per tile and month,
  so it should be read once per tile and month and cached.
- Advisory flags: 512 px blocks, LZW, no overviews.

## Maximum flood extent product

- It exists (PUM 3.1.1, released June 2023): a composite of the Observed Flood Extent
  layers of a period of **at most two months** for an area of interest.
- It is **generated on demand, not stored**: ordered in the portal or through
  `GET /download/max_flood_extent/{aoi_id}/{start_date}/{end_date}` and fetched with
  `GET /download/max_flood_extent/{order_id}`; the result arrives by e-mail link, as COG,
  GeoJSON or both. It needs a registered account and an area of interest saved in it.
- It is not an asset of the STAC items, so it cannot be read anonymously.
- For the comparison in step 3, a registered account is needed, or a file ordered by hand
  in the portal for the pilot area.

## Observation only

- The layers above are observations of Sentinel-1. The PUM lists topography (DEM, HAND) as
  auxiliary input of the algorithms and of the masks (exclusion of strong topography,
  correction of the reference water). The role inside the three algorithms is described
  in the Product Definition Document, which was not read here.
- No terrain-extended extent or depth is used: "Affected population" and "Affected land
  cover" and any third-party archive that extends GFM with a DEM stay out.

## Pilot event (September 2024, southern Poland)

- Search box 16.0–18.5° E, 50.0–51.2° N, 2024-09-10 to 2024-10-10: **116 items**,
  3–8 per acquisition day, on 20 days from 09-13 to 10-09 (about every 1–2 days).
- Scenes carry `anomaly_detected` and `flooded_pixels`, which allow screening before any
  file is opened.

## Consequences for the package

- The empty keys of `[sources]` can be filled with the verified values:
  `gfm_stac_url = "https://stac.eodc.eu/api/v1"`, `gfm_collection = "GFM"`,
  `gfm_assets = { flood = "ensemble_flood_extent", exclusion = "exclusion_mask",
  reference_water = "reference_water_mask", likelihood = "ensemble_likelihood" }`,
  `gfm_token_env = ""`.
- `gfm_reference_water_values` is `[1]` now; the mask also holds 2. Both should count as
  permanent or seasonal water: `[1, 2]`.
- The advisory flags are not used by the package yet; they mark scenes with unreliable
  contrast (frost, dry soil, wind) and are a candidate for excluding observations.

## Not verified

- The exact meaning of the reference water values 1 and 2.
- The maximum flood extent order itself (no account in this environment).
- How the three algorithms use HAND internally (Product Definition Document).
