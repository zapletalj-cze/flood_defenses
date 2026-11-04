<!-- Author: Jakub Zapletal -->
<!-- Date: 2026-11-04 -->
<!-- Podrobné vysvětlení metody balíku flood_defenses a zdůvodnění jednotlivých kroků (česky). -->

# Vrstva Defenses: jak funguje a proč

Tento dokument podrobně vysvětluje, jak balík `flood_defenses` odvozuje území chráněná
hrázemi a jejich dobu opakování ochrany, a především **proč** je každý krok navržen právě
takto. Spouštění, konfigurační klíče a výstupní soubory podrobně popisuje anglické
`README.md`; zde jsou zmíněny jen tam, kde pomáhají pochopit metodu. Čísla v textu jsou
naměřená v tomto vývojovém prostředí na syntetických datech (oddíl 13).

---

## Obsah

1. [Účel, vstupy a zásady](#1-účel-vstupy-a-zásady)
2. [Hlavní myšlenka: ochrana jako rozdíl dvou zaplavení](#2-hlavní-myšlenka-ochrana-jako-rozdíl-dvou-zaplavení)
3. [Bariéry: od linií hrází k rastru](#3-bariéry-od-linií-hrází-k-rastru)
4. [Terén bez hrází a HAND](#4-terén-bez-hrází-a-hand)
5. [Mezery v DTM: okraj pokrytí, nebo vodní hladina?](#5-mezery-v-dtm-okraj-pokrytí-nebo-vodní-hladina)
6. [Kompartmenty, vstupní buňka, zdroj a příznaky](#6-kompartmenty-vstupní-buňka-zdroj-a-příznaky)
7. [Velké DTM: koridor, jednotky, paralelismus, slučování](#7-velké-dtm-koridor-jednotky-paralelismus-slučování)
8. [Vodní zdroje a čerpací stanice z OpenStreetMap](#8-vodní-zdroje-a-čerpací-stanice-z-openstreetmap)
9. [Modelový odhad doby opakování](#9-modelový-odhad-doby-opakování)
10. [Pozorování I: satelitní záplavy GFM](#10-pozorování-i-satelitní-záplavy-gfm)
11. [Pozorování II: průtoky GloFAS a meze doby opakování](#11-pozorování-ii-průtoky-glofas-a-meze-doby-opakování)
12. [Validace](#12-validace)
13. [Výkon a škálování](#13-výkon-a-škálování)
14. [Co dokazují testy](#14-co-dokazují-testy)
15. [Omezení a otevřené otázky](#15-omezení-a-otevřené-otázky)
16. [Slovníček](#16-slovníček)

---

## 1. Účel, vstupy a zásady

### 1.1 Co vrstva je

Výsledkem jsou **kompartmenty**: polygony území chráněných ochrannými hrázemi. Každý
kompartment nese **dobu opakování ochrany** (RP, *return period*) a její interval
(`rp_low`…`rp_high`). Vrstva slouží k ořezu fluviálních hazardových map v
postprocessingu. Uvnitř kompartmentu s ochranou RP se zaplavení při povodních častějších
než RP neuvažuje.

### 1.2 Povolené a zakázané vstupy

| Povoleno | Zakázáno |
|---|---|
| LiDAR DTM (velké jednotlivé soubory 150–600 GB) | hazardové mapy |
| linie hrází (vlastní ML detekce a OSM) | výstupy hydraulických modelů |
| říční síť s plochou povodí | lokální hazardová data |
| volná data: OpenStreetMap, GHSL, Copernicus GFM, GloFAS | FABDEM |

**Proč bez hazardových map.** Vrstva se na hazard aplikuje. Kdyby z něj byla odvozena,
byla by kruhová: chyby hazardové mapy by se do ořezu přenesly dvakrát a vrstva by
nepřinesla žádnou nezávislou informaci. Proto stojí výhradně na terénu, hrázích, říční
síti a pozorováních. Pozorování (GFM, GloFAS) navíc vstupují jen jako *meze* doby
opakování; žádný kompartment nevytvářejí, neposouvají ani nemažou.

### 1.3 Priorita chyb: úplnost (recall)

Když je situace nejednoznačná, kompartment se **ponechá a označí příznakem**, místo aby se
zahodil.

- **Proč.** Vynechaný kompartment se v ořezané hazardové mapě projeví jako zaplavení
  chráněného území už při častých povodních. Za hrázemi typicky leží zástavba, takže je
  to chyba s velkým dopadem na škody.
- **Jak se ověřuje opačná chyba.** Falešný kompartment lze odhalit dalšími kroky (malá
  výška ochrany a z ní nízký strop RP, test „flooded“ ze satelitu, validace).

Všechna rozhodnutí, která vyměňují přesnost za úplnost, jsou v README výslovně označena.

### 1.4 Technické zásady a jejich důvody

- **Žádné parametry příkazové řádky, jen proměnné na začátku skriptu.** Každý běh je
  popsán souborem a několika proměnnými, které lze přepsat proměnnými prostředí
  `DEFENSES_<JMÉNO>` (používají je testy a dávkové běhy).
- **Všechny prahy, cesty, URL a CRS v jednom TOML souboru.** Výsledek je reprodukovatelný
  a v kódu nejsou skryté konstanty.
- **Kontrola konfigurace před výpočtem.** Každý skript nejdřív zkontroluje přesně klíče,
  které používá: typy, existenci souborů a vzájemnou konzistenci. Pak vypíše všechny
  problémy najednou a skončí. Chyba se tak neprojeví až po hodinách výpočtu.
- **Zástupné hodnoty (`placeholders`)** se při každém použití zapíší do logu jako
  varování. Nepotvrzená hodnota tak nemůže projít bez povšimnutí.
- **GDAL pro veškerý rastrový vstup a výstup** (bez rasteria a knihoven na něm závislých),
  **každý vektorový výstup jako samostatná jednovrstvá GeoPackage**. Obojí podle zadání;
  jednovrstvé soubory navíc vylučují záměnu vrstev.

---

## 2. Hlavní myšlenka: ochrana jako rozdíl dvou zaplavení

### 2.1 Princip jednou větou

> Hráz chrání území tehdy, když voda z řeky musí **s hrází vystoupat výš**, aby se do
> území dostala, než by musela bez ní.

Obě situace se počítají na téže ploše „výšky nad řekou“, takže se dají přímo odečíst.

### 2.2 HAND: výška nad nejbližší drenáží

`hand_undef` (*Height Above Nearest Drainage*) je výška buňky nad říční buňkou, do které
z ní stéká voda: `hand_undef = DTM_bez_hrází − river_ref`. Zde `river_ref` je výška říční
buňky, do které vede odtoková dráha (oddíl 4).

**Proč HAND, a ne absolutní výška.** Řeka má sklon. Hladina povodně je v prvním přiblížení
rovnoběžná s korytem, ne vodorovná. HAND převádí terén do souřadnic „výška nad místní
řekou“, takže jedna hodnota (např. 3 m nad řekou) znamená podobný povodňový stav po celé
délce toku. V absolutních výškách by území po proudu vždy vypadalo „níž“ a srovnání by
nedávalo smysl.

### 2.3 Úroveň přelití a priority flood

**Úroveň přelití** (*spill level*) buňky je nejnižší hladina nad řekou, při které voda z
libovolné říční buňky do buňky dosáhne. Voda se přitom pohybuje mezi čtyřmi sousedy:

```
spill(n) = min přes všechny cesty z řeky do n ( max výšky povrchu podél cesty ),   spill(řeka) = 0
```

Voda se do buňky dostane, jakmile hladina překročí nejvyšší bod (sedlo) na nejvýhodnější
cestě. Jde o úlohu minimaxové cesty. Řeší ji *priority flood*, tedy Dijkstrův algoritmus s
operátorem maxima: buňky opouštějí frontu v rostoucím pořadí úrovně přelití, takže
hodnota každé buňky je konečná hned při prvním vyjmutí. Implementace v `gis.py` (numba)
používá čtyřárnou haldu v poli a pro buňky na aktuální úrovni nebo pod ní obyčejnou frontu
(vylepšený priority flood podle Barnese). Na jednotce 30 M buněk trvá jedno zaplavení
11 s.

### 2.4 Dva povrchy, dvě zaplavení, jeden rozdíl

| Veličina | Povrch | Význam |
|---|---|---|
| `spill_undef` | `hand_undef` | kdy by se voda do buňky dostala bez hrází |
| `spill_def` | `hand_undef`, ale na buňkách bariér `max(koruna − river_ref, hand_undef)` | kdy se tam dostane s hrázemi |
| `margin` | `spill_def − spill_undef ≥ 0` | o kolik hráze zvyšují potřebnou hladinu |

Buňka je **chráněná**, pokud platí obě podmínky:

- `margin ≥ min_margin_m` (výchozí 0,3 m). Rezerva pár centimetrů není skutečná ochrana,
  ale šum DTM a nejistota výšky koruny.
- `spill_undef ≤ max_level_m`. Území dosažitelné až při extrémní hladině nad řekou není
  fluviální inundace; zaplavené by nebylo ani bez hrází.

Kompartmenty jsou souvislé (4-sousedství) oblasti chráněných buněk o ploše aspoň
`min_area_m2`. **Úroveň ochrany** `protection_level_m` je nejmenší `spill_def` v
kompartmentu: výška hladiny nad řekou, při které voda do kompartmentu poprvé vnikne
(přeteče přes nejnižší místo hráze na nejvýhodnější cestě). Tato výška později určuje
hloubky vody i strop RP (oddíl 9).

**Příklad.** Koruna hráze leží 4 m nad řekou, terén za ní 2 m nad řekou.
- Bez hráze: `spill_undef = 2 m`.
- S hrází: `spill_def = 4 m`.
- Rezerva je tedy 2 m, buňka je chráněná a `protection_level_m = 4 m`.

### 2.5 Důsledky, které na první pohled překvapí, ale jsou správné

- **Části okruhu výš než jeho nejnižší koruna nejsou chráněné.** Voda by k nim dosáhla
  při své vlastní výšce tak jako tak. Hráz je tedy chrání jen tam, kde opravdu mění, kdy
  se voda dostane dovnitř.
- **Vnořené okruhy.** Vnitřní kompartment má jako úroveň ochrany vyšší z obou korun. Test
  to ověřuje v obou pořadích výšek.
- **4-sousedství.** Voda nesmí projít mezi dvěma buňkami hráze, které se dotýkají jen
  rohem. Rastr se proto navíc „utěsňuje“ (oddíl 3).

---

## 3. Bariéry: od linií hrází k rastru

1. **Sloučení zdrojů.** Všechny soubory hrází se načtou, převedou do CRS DTM a rozdělí na
   jednoduché linie. Atribut `source` nese původ (ML, OSM).
   - **Kdy se to dělá.** Jednou za zemi, ne v každé jednotce: všechny jednotky pak vidí
     identické bariéry a výsledky na švech sedí.
2. **Přemostění mezer.**
   - **Problém.** ML detekce i OSM přerušují linie u mostů, cest a propustků. Skutečná hráz
     je tam obvykle spojitá (nájezd, mobilní hrazení).
   - **Pravidlo.** Konec linie, u kterého je do `gap_tolerance_m` (20 m) jiná linie, se
     spojí: s jejím koncem, pokud je v toleranci, jinak s nejbližším bodem (T-spoj).
   - **Výjimky.** Konce blíž než `touch_tolerance_m` se považují za spojené. Dlouhá linie,
     jejíž vlastní konce jsou v toleranci, se uzavře do okruhu.
   - **Proč.** Úplnost. Skutečně otevřená mezera širší než tolerance zůstane otevřená a
     území za ní chráněné není; to ověřuje test se širokou mezerou.
3. **Rasterizace s `ALL_TOUCHED`**, aby tenká linie neměla v rastru díry. Pak se
   **uzavřou čistě diagonální kontakty**: přidá se vyšší ze dvou rohových buněk.
   Řeka křížící „schodovitou“ hráz tak vždy přetne buňku bariéry.
4. **Výška koruny z nativního DTM.**
   - **Postup.** Bere se maximum nativních pixelů v okruhu `crest_search_radius_m` od
     linie. Čte se jen v blocích kolem hrází.
   - **Proč nativní rozlišení.** Průměrování na pracovní rozlišení 5 m korunu snižuje.
   - **Proč maximum v okolí.** Geometrie linie bývá od skutečné koruny posunutá o metry.
   - **Proč bez overviews.** Čte se s `resampleAlg="max"` a **vypnutými overviews**,
     protože overviews počítané průměrem by korunu snížily.
   - **Efektivní poloměr** je až o jednu buňku větší než nastavený; to jde ve prospěch
     úplnosti.
5. **Přemostěné buňky** dostanou korunu lineárně interpolovanou mezi korunami na koncích
   mostu.
6. **Řeka křížící hráz.** Ve výchozím stavu „vyhrává řeka“: buňky křížení nejsou
   bariérou, protože řeka linií hráze protéká. Přítoky uzavřené stavidlem nebo
   přečerpávané řeší krok OSM (oddíl 8).

---

## 4. Terén bez hrází a HAND

### 4.1 DTM bez hrází

Buňky bariér rozšířené o `barrier_removal_width_m` se nastaví na nodata a vyplní
interpolací `gdal.FillNodata`.

- **Proč.** HAND musí popisovat terén *bez* hráze, jinak by se těleso hráze tvářilo jako
  vyvýšený terén. Účinek hráze vstupuje do výpočtu jen přes druhý povrch (`spill_def`).
- **Proč pás, a ne jen linie.** Těleso hráze je širší než jeho osa.
- **Oprava B2.** Říční buňky v tomto pásu si ponechají svou výšku, protože jsou korytem.
  Interpolují se jen tam, kde řeka hráz kříží a DTM ukazuje těleso hráze. Původní verze
  přepisovala dno koryta interpolovanými výškami břehů.

### 4.2 Výpočet HAND

1. **Říční buňky.** Řeky, které jsou zdrojem vody, se rasterizují po ose (bez
   `ALL_TOUCHED`), aby přetínaly co nejméně buněk bariér.
2. **Směrování odtoku** běží na **upravené kopii** DTM bez hrází:
   - **Zahloubení koryta s náběhem (oprava B1).** Říční buňky se sníží o
     `stream_burn_m`, lineárně ubývající k nule ve vzdálenosti `stream_burn_buffer_m`.
     Bez toho vedly odtokové dráhy D8 *podél* koryta až k jeho nejnižšímu bodu na kraji
     okna. Celé okno pak dostalo jako referenci výšku řeky z nejníže položeného místa a
     HAND vyšel nadhodnocený o sklon × délku okna.
   - **Úplné prolomení (breaching).** Priority flood z říčních buněk: bezodtoká deprese se
     „prokope“ směrem k výtoku. Ploché oblasti dostanou nepatrný spád (1e-5 m).
   - **D8 nejstrmějším spádem.** Na upraveném povrchu je příjemce vždy níž, cesty jsou
     acyklické a říční buňka se pro všechny buňky určí v O(n).
3. **`river_ref`** je výška dosažené říční buňky z **neupraveného** DTM bez hrází.
   Úprava slouží jen k nalezení cesty, ne k výškám.
4. **Záloha.** Buňky bez spojení s řekou dostanou nejbližší říční buňku; jejich počet se
   zapíše do logu.
5. **Výtoky jsou jen modelované řeky, které jsou zdrojem vody.** Boční povodí se tak
   směrují k řece, kterou hazardové mapy obsahují. Každá buňka má referenci vůči řece,
   z níž by k ní povodeň přišla.

---

## 5. Mezery v DTM: okraj pokrytí, nebo vodní hladina?

LiDAR mívá nodata na vodních hladinách. Nodata však vzniká ze dvou velmi odlišných
důvodů:

| Typ | Co to je | Zacházení | Proč |
|---|---|---|---|
| chybějící pokrytí | území mimo snímkování | zůstane nodata a funguje jako zeď; kompartment u něj se ponechá s příznakem `touches_dtm_edge` (oprava B3) | terén je neznámý; recall |
| mezera uvnitř pokrytí | vodní hladina, výpadek dat | interpoluje se | řeka zapsaná jako nodata by neměla říční buňky, tedy ani zdroj vody, ani kompartmenty |

**Jak je rozlišit.** Pokrytí je **morfologické uzavření** platných buněk s poloměrem
`max_fill_distance_m`, spočtené dvěma průchody vzdálenostní transformací. Cena proto
nezávisí na poloměru.
- **Co uzavření dělá.** Zavře mezery užší než dvojnásobek poloměru, včetně místa, kde řeka
  vytéká z pokrytí, ale přímý okraj pokrytí nepohne.
- **Okraj okna.** Maska se před uzavřením obloží „žádnými daty“, aby se okraj pokrytí
  blízko hranice okna neuzavřel proti ní.

**Vodní hladiny (`cap_river_holes`)**, tedy mezery, které křižuje řeka:
- **Říční buňky** se interpolují **podél osy řeky** mezi nejbližšími platnými říčními
  buňkami (podélný profil).
- **Ostatní buňky takových mezer** a případ, kdy je celá řeka v okně nodata: výplň se
  shora omezí nejnižší platnou buňkou v dosahu, tedy úrovní hladiny na okraji.
- **Proč.** Interpolace jen z břehů leží *nad* hladinou. V testu zvedla `river_ref` asi o
  1 m a vytvořila falešný kompartment. Test T3 udělá z celé řeky nodata a ověřuje, že
  výsledek je stejný jako s úplným DTM.

---

## 6. Kompartmenty, vstupní buňka, zdroj a příznaky

### 6.1 Vstupní buňka a zdroj

Priority flood vrací pro každou buňku i **výchozí říční buňku** (odkud vede její
minimální cesta) a **pořadí dosažení**.

- **Vstupní buňka (`entry_x`, `entry_y`):** buňka kompartmentu, které zaplavení s hrázemi
  dosáhlo jako první, tedy místo, kudy voda vnikne.
  - **Proč pořadí, a ne minimum `spill_def`.** Všechny buňky pod korunou mají stejnou
    úroveň přelití, takže minimum je nejednoznačné; pořadí dosažení to rozhodne.
- **Zdroj (`source_x`, `source_y`):** výchozí říční buňka vstupní buňky, tedy řeka, jejíž
  voda vnikne první. Podle ní se určuje plocha povodí (oddíl 9) a hledá buňka GloFAS
  (oddíl 11).

### 6.2 Příznaky místo mazání

| Příznak | Význam | Proč příznak |
|---|---|---|
| `truncated` | dotýká se okraje okna jednotky, může pokračovat dál | řeší se při slučování (7.4) |
| `touches_corridor_edge` | sousedí s buňkami mimo zpracovávanou doménu | výsledek může být neúplný |
| `touches_dtm_edge` | sousedí s chybějícím pokrytím, které ho uzavírá jako zeď | ochrana může být umělá; recall |
| `has_pumping` | uvnitř nebo do `levee_adjacent_distance_m` je čerpací stanice | informace pro odhad RP a kontrolu |
| `no_levee` | u hranice není žádná linie hráze | území uzavřené terénem; RP jen ze stropu |

---

## 7. Velké DTM: koridor, jednotky, paralelismus, slučování

DTM o velikosti stovek GB se nikdy nečte celé. Všechna čtení jdou přes okenní funkce
`gis.Raster` na výpočetních jednotkách.

### 7.1 Stopa a koridor

- **Stopa platných dat** se čte z DTM decimovaného na `footprint_resolution_m` jedním
  voláním. GDAL použije overview, pokud existuje, jinak dekóduje jen pásy obsahující
  vzorkované řádky. Na testovacím DTM to trvá 0,1 s na soubor.
- **Koridor** jsou řeky obalené `corridor_buffer_m`, oříznuté stopou.
  - **Proč.** Území daleko od modelovaných řek jimi zaplaveno být nemůže. Práce se tak
    omezí na podíl údolních niv, typicky zlomek plochy státu.

### 7.2 Uložení DTM: proč na něm záleží

U pásového uložení (výška bloku pod 64 řádků nebo blok přes celou šířku) bez overviews
dekóduje každé okenní čtení všechny nativní pixely. Naměřené časy čtení terénu jedné
jednotky 30 M buněk:

| Uložení | Čtení terénu |
|---|---|
| pásové, bez overviews | 8–9 s |
| po konverzi A (dlaždice 512², ZSTD, overviews) | 1,7–2,3 s |

- **`prepare_dtm.py` nic nepředpokládá, ale měří.** Hlásí blokovou strukturu, kompresi a
  overviews a spustí benchmark: 20 náhodných oken 4096², stejná okna přes `gdal.Warp` a
  sekvenční čtení.
- **Použití overview se ověřuje.** Z ladicího výpisu GDAL se čte zdrojové okno; je-li
  menší než nativní, overview bylo použito.
- **Konverze** (A celý soubor, B jen koridor s lemem) se spustí jen na vyžádání. Předem
  vypíše odhad místa na disku a času z naměřených hodnot a zkontroluje volné místo.

### 7.3 Jednotky: jádro, lem, doména

- **Jádro.** Čtverec `unit_core_size_m` zarovnaný na pracovní mřížku.
- **Okno.** Jádro rozšířené o lem `unit_halo_m`.
- **Doména.** Uvnitř okna se počítají jen buňky koridoru s lemem; ostatní jsou nodata a
  všechny algoritmy je přeskočí.

**Pravidlo `unit_halo_m ≥ corridor_buffer_m` a jeho zdůvodnění:**
1. Buňka jádra může ležet až `corridor_buffer_m` od své řeky.
2. Okno přesahuje jádro o lem na všechny strany.
3. Nejbližší bod řeky tak je pro každou buňku jádra uvnitř okna jen tehdy, když lem ≥
   obalová vzdálenost.
4. Jinak by HAND neměl platnou referenci. Konfigurace proto kratší lem odmítne.

**Paměťový model.** RAM jednotky ≈ buňky okna × bajty na buňku.
- **Naměřeno** pro `compartments.py`: 93–95 B na buňku.
- **Kontrola.** `build_units.py` se zastaví, pokud by největší okno překročilo
  `max_ram_gb_per_worker`, a vypíše největší jádro, které se vejde.

### 7.4 Paralelní běh a kontrolní body

**Fond procesů** se spouští metodou `spawn` a s jedním úkolem na proces.
- **Proč nový proces pro každou jednotku.**
  - Špička paměti se měří přesně po jednotkách.
  - Paměť se po každé jednotce vrátí systému (numpy a GDAL ji jinak drží).
  - Linux a Windows se chovají stejně.
- **Kernely numba** se přeloží jednou v hlavním procesu a pracovníci je načtou z cache.
- **Pořadí.** Jednotky jdou od největší, odhad zbývajícího času se počítá z hotových
  buněk.
- **Kontrolní body.** Každá jednotka končí souborem `done.json` (časy, špička paměti,
  počty) zapsaným atomicky. Opakovaný běh přeskočí hotové jednotky. Pád jednoho
  pracovníka, např. pro nedostatek paměti, se ohlásí a opakuje se jen ta jednotka.

### 7.5 Slučování jednotek

Okna se překrývají o dvojnásobek lemu, takže stejný kompartment může najít více jednotek.

1. **Skupiny.** Kandidáti různých jednotek, kteří se překrývají aspoň
   `merge_overlap_ratio` menšího z nich, tvoří skupinu (souvislé komponenty grafu).
2. **Skupina s úplným kandidátem.** Ponechá se úplný (neuseknutý) kandidát, jehož
   reprezentativní bod leží v jádru vlastní jednotky, při shodě ten nejdál od okraje okna.
   Okruh přes hranici jader se tak ponechá právě jednou (test švu).
3. **Skupina useknutá ve všech oknech** (kompartment delší než okno):
   - **Geometrie.** Části se **sjednotí** (recall).
   - **Úroveň ochrany, vstup a zdroj** se převezmou z části s **nejnižší** úrovní ochrany.
   - **Proč nejnižší.** Okraj okna funguje jako zeď, takže ochranu může jen nadhodnotit,
     nikdy podhodnotit.
4. **Pravidlo rozporu.**
   - **Problém.** Okraje oken mohou vytvořit umělou ochranu, například území za dlouhou
     hrází, které vypadá uzavřené jen díky okrajům.
   - **Pravidlo.** Kompartment useknutý ve všech oknech se zahodí, pokud jednotky, jejichž
     okna pokrývají aspoň `contradiction_ratio` jeho plochy, na jeho místě nenašly nic.
     Tyto jednotky vidí z hráze víc a „mlčí“.
5. **Identifikátory.** `compartment_id` = kód země + reprezentativní bod zaokrouhlený na
   10 m, takže je stabilní mezi běhy.

---

## 8. Vodní zdroje a čerpací stanice z OpenStreetMap

### 8.1 Problém

Přítok procházející okruhem hrází je obvykle za povodně **uzavřený** stavidlem v místě
křížení s hrází, nebo se voda z vnitřní strouhy **přečerpává** přes hráz. Hydraulicky
pak není zdrojem vody pro chráněné území. Bez této informace „vyhrává řeka“: přítok se
považuje za otevřený, voda jím vtéká dovnitř a okruh vyjde **nechráněný**. To je ztráta
úplnosti právě u nejtypičtějších polderů.

### 8.2 Řešení

1. **`download_osm.py`** stáhne výřez OSM země. Jedním SQL dotazem přes body, linie a
   plochy (vlastní `osmconf.ini`) vytáhne hydraulické objekty (stavidla, vrata, jezy,
   čerpací stanice; plochy jako body uvnitř) a odvodňovací kanály. Zapíše počty podle
   hodnot tagů, aby šly seznamy tagů ladit.
2. **`osm_hydraulics.py`** označí říční prvek jako **nezdroj** (`is_source = False` s
   důvodem v `reason`), pokud:
   - u některého jeho křížení s hrází leží do `gap_tolerance_m` uzavírací objekt
     (`osm_tags_closing`, výchozí stavidla a vrata), nebo
   - u jeho **dolního konce** leží čerpací stanice.
3. **Dolní konec z výšek DTM.**
   - **Postup.** Medián nativního DTM v okolí obou konců linie; nižší konec je dolní.
   - **Proč z DTM.** Směr digitalizace říčních linií není spolehlivý; test obsahuje linii
     digitalizovanou proti proudu.
   - **Kdy rozhodne směr linie.** Jen při rovnosti nebo chybějících výškách; sloupec
     `downstream_method` uvádí, co rozhodlo.
4. **Uzavřené toky se vyřadí** ze zdrojů vody **i ze sítě pro HAND**.
   - **Proč i z HAND.** Jinak by HAND uvnitř kompartmentu měřil výšku nad vnitřní strouhou
     místo nad hlavní řekou a hloubky i úroveň ochrany by byly chybné.
   - **Důsledek v testu.** Okruh s přítokem uzavřeným stavidlem je chráněný až po kroku OSM.
5. **`osm_pumping.gpkg`** obsahuje čerpací stanice u hrází nebo řek. Z nich se odvozuje
   příznak `has_pumping`.

**Omezení.** Rozhoduje se o celém prvku, jak je v datech. Prvek pokračující daleko za
hrází se uzavře celý, proto je vhodné mít říční síť rozdělenou na soutocích.

---

## 9. Modelový odhad doby opakování

### 9.1 Proč model

Návrhovou úroveň ochrany nelze spočítat hydraulicky, protože hazardová data jsou
zakázaná. Odhaduje se proto z toho, **co hráz chrání** (expozice) a **jak velká je řeka**.
Pořadí kompartmentů podle těchto veličin se pak promítne na předpokládané národní
rozdělení tříd ochrany.

### 9.2 Hodnota expozice (GHSL)

1. **Objemy zástavby.** Z GHSL (celkový a nebytový objem budov, 100 m, Mollweide) se
   spočte obytný objem = celkový − nebytový, oříznutý nulou. Podlahová plocha = objem /
   výška podlaží.
2. **Hodnota** = podlahová plocha × maximální škoda na m² podle země (tabulka
   `max_damage_csv` z JRC, dodaná vámi).
3. **Převzorkování se zachováním součtů.** Hodnota na m² se převzorkuje průměrem a
   vynásobí plochou pracovní buňky.
   - **Proč to funguje.** Mollweide je stejnoplochá projekce, takže součet zůstává
     zachován.
   - **Kontrola.** V každé jednotce se součet porovná s nezávislým součtem, ve kterém jsou
     zdrojové buňky váženy **přesným** podílem plochy uvnitř okna. Naměřená shoda je
     0,03–0,05 % (tolerance 1 %).
4. **Rozsah výpočtu.** Hodnota se počítá jen v obálce kompartmentů jednotky: 2–5 s místo
   15 s.

### 9.3 Škoda při selhání ochrany

- **Hloubka v buňce** = `max(0, protection_level_m − hand_undef)`.
  - **Proč.** Je to hladina v okamžiku přelití vyjádřená v HAND, tedy hloubka, kterou by
    území mělo, kdyby voda vystoupala na úroveň ochrany a hráz byla přelita.
- **Podíl škody** z křivek hloubka–škoda (`depth_damage_csv`, obytné a nebytové),
  lineárně interpolovaný.
- **V** = Σ hodnota × podíl škody přes kompartment: škoda, které hráz při své návrhové
  události zabraňuje.

### 9.4 Skóre a třídy

- **L** = délka linií hrází u hranice kompartmentu. Kompartment bez hrází (`no_levee`)
  skóre nedostane a RP má jen ze stropu.
- **A** = plocha povodí u zdrojové říční buňky (z atributu řeky nebo z rastru).
- **Skóre** = `w_v × pořadí(log V/L) + w_a × pořadí(log A)`, pořadí jako percentil v rámci
  země.
- **Proč V/L.**
  - **Logika přínosů a nákladů.** Přínos ochrany je zabráněná škoda V, náklady rostou s
    délkou hráze L. Vysoká hodnota na metr hráze ospravedlňuje vyšší standard.
  - **Proč pořadí.** Nevyžadují kalibrované absolutní hodnoty a jsou odolná vůči
    extrémům i jednotkám.
- **Proč A.** Velké řeky mívají vyšší návrhové standardy: státní správa, větší
  investice, delší historie ochrany.
- **Kvantilové mapování.** Kompartmenty seřazené podle skóre dostanou třídy
  `distribution_classes` tak, aby podíly odpovídaly `distribution_shares`; stejné skóre
  dá stejnou třídu.
  - **Proč.** Kolik hrází má jaký standard, se dá odhadnout národně (zástupná hodnota k
    potvrzení), ale ne, které to jsou. Model dodává právě pořadí.

### 9.5 Strop, interval, ruční zásahy

- **Strop `rp_cap`** je schodovitá funkce úrovně ochrany: první hodnota `cap_levels_m`
  **nad** úrovní vybere odpovídající `cap_rp`.
  - **Proč.** Fyzikální věrohodnost. Hráz vystupující 0,3 m nad údolní nivu neudrží
    stoletou povodeň.
- **Výsledek.** `rp = min(rp_model, rp_cap, max_rp)`. `rp_source` je `model`, `cap` nebo
  `override`.
- **Interval.** `rp_low` je předchozí třída a `rp_high` následující třída omezená stropem
  a `max_rp`, nikdy pod `rp`.
  - **Proč ±1 třída.** Vyjadřuje nejistotu modelu poctivěji než bodová hodnota a dává
    pozorováním (oddíl 11) prostor interval zúžit.
- **Ruční zásahy (overrides).** Bod se známou RP uvnitř kompartmentu má vždy přednost:
  `rp = rp_low = rp_high`. Pozorování ho nikdy nemění.

### 9.6 Doplňková expozice z land cover (modul HAND_YZ)

**Proč.** GHSL měří budovy. Hodnota mimo budovy proto v 9.2 chybí:
- plodiny;
- silnice, železnice a přístavy;
- venkovní zařízení průmyslových areálů.

Modul `HAND_YZ.py` ji doplní z dat o využití a krytu území. Hodnotu z GHSL tedy
**doplňuje**, nenahrazuje.

1. **Data.**
   - **ESA WorldCover 10 m** (2021, v200) všude. `download_worldcover.py` stáhne dlaždice
     oficiální mřížky 3° × 3°, které zasahují do koridoru, a sestaví `worldcover.vrt`.
   - **Urban Atlas** tam, kde ho dodáte (ruční stažení). Má přednost všude, kde leží jeho
     polygony, i s třídami bez hodnoty, protože je podrobnější a rozlišuje využití
     území.
   - **Vícevrstvý soubor** se převede na jednovrstvou kopii přes `gdal.VectorTranslate`.
2. **Kategorie.** Třídy WorldCover a Urban Atlasu se přiřadí kategoriím expozice
   (např. `agriculture`, `built_up`, `transport`). Každá kategorie má:
   - hodnotu na m² pozemku (`lulc.max_damage_csv`);
   - vlastní křivku hloubka–škoda (`lulc.depth_damage_csv`).
3. **Škoda.** Hloubka je stejná jako v 9.3: `max(0, úroveň ochrany − HAND)`.
   Hodnota buňky = plocha × hodnota na m² její kategorie; škoda = hodnota × podíl škody
   z křivky kategorie. Součty za kompartment a kategorii se zapíší do
   `exposure_lulc.csv`.
4. **Spojení.** S `lulc.enabled = true` počítá `rp_index.py`
   `value_exposed = value_ghsl + value_lulc`. Skóre pak vychází ze součtu a obě části
   zůstanou ve výstupu.
5. **Pozor na dvojí započtení.** Hodnoty na m² smějí pokrýt jen to, co nejsou budovy z
   GHSL.
   - **Třída „built-up“ WorldCoveru** obsahuje i budovy, proto má mít hodnotu jen za
     komunikace, sítě a venkovní majetek.
   - **Nejisté hodnoty** nastavte na 0.
6. **Kontrolní body.** Každá jednotka si uloží otisk (hash) vstupů: hodnot, křivek,
   převodníků tříd a dat souborů. Opakovaný běh přepočítá právě jednotky, kterým se
   vstupy změnily (testováno).

**Co je ověřené.**
- **WorldCover:** vzor URL a mřížka z oficiálního stahovacího skriptu ESA, skutečný
  požadavek na dlaždici, kódy tříd a licence CC BY 4.0 přímo z metadat dlaždice.
- **Urban Atlas:** název pole a seznam kódů ověřit nešlo (hostitelé EEA byli v této
  session blokovaní), proto jsou v konfiguraci prázdné.

---

## 10. Pozorování I: satelitní záplavy GFM

### 10.1 Co je GFM

Copernicus Global Flood Monitoring vyhodnocuje každý snímek radaru Sentinel-1 (SAR, od
roku 2015) souborem tří algoritmů. Vrstvy produktu:
- **rozsah záplavy:** 0/1 v rozlišení 20 m;
- **maska vyloučení:** místa, kde radar vodu spolehlivě nevidí, např. hustá vegetace,
  zástavba, radarový stín;
- **referenční (trvalá) voda.**

Názvy vrstev v API a jejich kódy pixelů zatím nejsou ověřené (oddíl 15).

### 10.2 Stažení (`download_gfm.py`): proč takhle

1. **Jedno vyhledávání STAC pro všechny jednotky** s průchodem stránkami (GET i POST
   odkazy `next`) a volitelným tokenem z proměnné prostředí. Token se nikdy neukládá do
   konfigurace ani do logu.
2. **Předběžná kontrola: jedno čtení na scénu.** Vrstva záplavy se jednou přečte v
   rozlišení 100 m (průměr = podíl zaplavení) přes okna všech jednotek, které scéna
   pokrývá. U každé jednotky se sečte zaplavená plocha v její doméně.
   - **Proč jedno čtení na scénu.** Země má stovky jednotek a každá jednotka asi tisíc
     snímků Sentinelu-1 od roku 2015. Čtení zvlášť pro každou jednotku by znamenalo
     statisíce síťových požadavků. Takto jich je o tolik méně, kolik jednotek scéna
     pokrývá.
   - **Proč je to levné.** COG obslouží čtení ve 100 m z overviews, tedy kilobajty.
   - **Pravidlo.** Okno jednotky s nejvýš `min_valid_pixels` zaplavenými buňkami nemůže
     dát platný test a přeskočí se. Většina snímků je suchá, takže se kopíruje jen málo
     scén.
3. **Kopie oken** vrstev záplavy, vyloučení a referenční vody (a volitelně
   pravděpodobnosti) zůstávají v mřížce a CRS GFM.
   - **Proč.** Převzorkuje se jen jednou, přímo na mřížku jednotky v dalším kroku;
     dvojí převzorkování by posunulo hrany 20m pixelů.
4. **Vlákna (`gfm_download_threads`).** Práce čeká na síť a GDAL během čekání uvolňuje
   GIL, takže stačí vlákna místo procesů.
5. **Kontrolní bod `screening.csv`.** Obsahuje každou prověřenou dvojici scéna–jednotka.
   Opakovaný běh, třeba za rok, stáhne jen nové snímky.

### 10.3 Test kompartmentu (`gfm_evidence.py`)

Suchý kompartment v suchý den nedokazuje nic. Každá scéna proto odpovídá na srovnávací
otázku: **byl kompartment suchý, zatímco nechráněná niva ve stejné výšce hned vedle byla
zaplavená?** Teprve srovnání ve stejné scéně dělá ze suchého kompartmentu důkaz.

1. **Platná buňka** má data, není vyloučená, není trvalá voda a má HAND. **Mokrá buňka**
   je platná a zaplavená.
2. **Uvnitř:** podíl mokrých platných buněk kompartmentu (`inside_wet_frac`).
3. **Srovnávací zóna (`comp_wet_frac`)** zahrnuje platné buňky, které:
   - leží mimo všechny kompartmenty,
   - jsou do `comparison_distance_m` od kompartmentu,
   - mají `hand_undef` nejvýš medián HAND uvnitř kompartmentu.

   Podmínka HAND omezí srovnání na terén, kam by nechráněná povodeň do úrovně
   kompartmentu také došla. Jsou-li tyto buňky mokré, voda u hráze stála aspoň tak
   vysoko.
4. **Výsledek testu:**
   - `no_test`, pokud chybí dost platných buněk uvnitř nebo v zóně, nebo pokud zóna není
     dost mokrá (`f_wet_comparison`);
   - `held` („vydržel“) při podílu uvnitř ≤ `f_dry_inside` (0,05: tolerance šumu SAR a
     loužek);
   - `flooded` („zaplaven“) při podílu ≥ `f_wet_inside` (0,3);
   - jinak `inconclusive`.
5. **Výstupy.**
   - `gfm_evidence.csv`: řádek pro každou dvojici kompartment × scéna.
   - `gfm_evidence.gpkg`: tytéž řádky s polygonem kompartmentu pro kontrolu v QGIS.
   - **Jedna povodeň.** S `DATE_FROM` a `DATE_TO` se zpracuje jen jedna událost; takto
     uděláme kontrolu, až určíte povodeň.

### 10.4 Upozornění na chybějící bariéry (`HINTS = True`)

1. **Hladina scény v HAND** je 90. percentil HAND mokrých buněk **spojených s řekou**.
   - **Spojitost** se počítá přes mokré buňky, trvalou vodu a osy řek.
   - **Proč i trvalá voda.** Řeka sama je trvalou vodou, takže ve vrstvě záplavy nikdy
     není „mokrá“. Bez toho by se záplava od své řeky oddělila; tuto chybu našel test.
2. **Testované buňky** splňují všechny podmínky:
   - jsou platné a mimo kompartmenty,
   - leží aspoň `hint_exclusion_m` od hrází a kompartmentů,
   - jejich HAND je pod hladinou − `hint_margin_m`.
3. **Upozornění** je shluk buněk, které zůstaly suché aspoň v `hint_dry_share` testovaných
   scén, při aspoň `hint_min_scenes` testech a ploše aspoň `min_area_m2`.
4. **Proč `hint_exclusion_m` (doplněk k zadání).**
   - **Problém.** Hráz, která vydržela, je sama suchá; suché jsou i její paty a smíšené
     20m pixely. Tyto buňky leží v `hand_undef` pod hladinou, protože hráz byla z povrchu
     odstraněna. Bez vyloučení by každá hráz, která vydržela, vypadala jako „chybějící
     bariéra“ podél svého tělesa.
   - **Test.** Bez vyloučení vzniknou 3 upozornění místo 1.
5. **Omezení.** HAND je relativní k místní drenáži. Při zpětném vzdutí velké řeky do
   údolí přítoku mají buňky u přítoku nízký HAND vůči přítoku, ale hladinu určuje hlavní
   řeka. Mohou zůstat legitimně suché a objevit se jako upozornění. Upozornění jsou proto
   **kandidáti ke kontrole**, ne výsledek, a žádný kompartment nemění.

---

## 11. Pozorování II: průtoky GloFAS a meze doby opakování

### 11.1 Myšlenka

Snímek SAR říká, **zda** byl kompartment mokrý. Neříká, **jak vzácná** byla povodeň. Tu
informaci dává průtok v řece ve stejný den, srovnaný s dlouhou řadou ročních maxim:

- **„vydržel“ (`held`)** při události s dobou opakování T → ochrana je **aspoň** T
  (dolní mez);
- **„zaplaven“ (`flooded`)** při události s dobou opakování T → ochrana je **nejvýš** T
  (horní mez).

### 11.2 Data (`download_glofas.py`)

- **Řada.** Historická reanalýza GloFAS (`cems-glofas-historical`, ověřeno na stránce
  datasetu): denní průtok `dis24`, mřížka 0,05°, od roku 1979.
- **Plocha povodí.** Statická mapa plochy povodí na stejné mřížce.
- **Dvě cesty.**
  - **Lokální soubory.** GRIB nebo NetCDF stažené formulářem EWDS.
  - **Roční stahování přes `cdsapi`.** Jen pokud ho schválíte a nainstalujete.
- **Čtení jen přes GDAL.** U GRIB je čas pásma `GRIB_VALID_TIME`. U NetCDF je to dimenze
  času s atributem `units` („days/hours/seconds since …“).
- **Přiřazení dne.** Denní průměr označený koncem svého 24h období patří ke dni, který
  končí. Přepínač `glofas_valid_time_ends_day` je zástupná hodnota k ověření; posun o den
  je důležitý hlavně u bystřin.
- **Výstup.** Řady buněk s dostatečnou plochou povodí v kompaktním `glofas_discharge.npz`.

### 11.3 Přiřazení buňky (`rp_bounds.py`)

Ke zdrojové buňce kompartmentu se hledají buňky GloFAS:
1. **Okruh.** Do `glofas_search_radius_km` od zdrojové buňky.
2. **Filtr.** Poměr ploch povodí v rozsahu `upstream_ratio_min`…`upstream_ratio_max`
   (0,7–1,3).
3. **Volba.** Vyhrává poměr nejbližší 1 (v logaritmu).

**Proč poměr ploch, a ne vzdálenost.** V mřížce 0,05° (asi 5 km) je nejbližší buňka často
na jiné řece. Poměr ploch vybere buňku, která reprezentuje **stejný úsek toku**. Test
obsahuje druhou buňku v okruhu s poměrem 0,92; vyhraje buňka s poměrem 1.

### 11.4 Frekvenční analýza

Roční maxima přiřazené buňky (jen roky s aspoň `glofas_min_days_per_year` dny) se
proloží Gumbelovým rozdělením metodou L-momentů:

```
b0 = průměr,  b1 = Σ (i−1)/(n−1) · x(i) / n  (vzestupně),  l1 = b0,  l2 = 2·b1 − b0
α = l2 / ln 2,   ξ = l1 − 0,5772 · α
RP(Q) = 1 / (1 − exp(−exp(−(Q − ξ)/α)))
```

- **Proč Gumbel a L-momenty.** Pro zhruba 40 ročních maxim je odhad stabilní a nemá
  parametr tvaru, který by se dal přeučit.
- **Kontrola správnosti.** Medián ročních maxim má RP ≈ 2, protože medián má
  pravděpodobnost nepřekročení 0,5. Test to splňuje (2,00 a 2,003).
- **Slabina: těžké chvosty.** Chvost Gumbelova rozdělení je lehčí než u řek s těžkým
  chvostem. Velký průtok tak dostane vyšší RP, než ve skutečnosti má.
  - **Dopad na dolní meze.** „Vydržel“ dává příliš vysoké dolní meze, což je u dolní meze
    nebezpečný směr.
  - **Dopad na horní meze.** „Zaplaven“ dává příliš vysoké horní meze, což je jen uvolní.
  - **Proč ne GEV.** GEV se třemi parametry by chvost zachytil, ale pro jednotlivou buňku
    se 40 maximy je nestabilní. Dalším krokem by bylo regionální sdružování buněk.

### 11.5 RP události a pozorované meze

- **`held`.** RP průtoku **v den snímku** je dolní mez: hráz vydržela aspoň tolik.
- **`flooded`.** RP **maxima** průtoku za `flooded_window_days` do dne snímku včetně je
  horní mez. Kompartment zaplavila událost nejvýš této velikosti a okno zachytí kulminaci,
  která mohla přejít před přeletem družice.
- **Za kompartment:**
  - `obs_rp_low` = maximum přes testy „vydržel“;
  - `obs_rp_high` = minimum přes testy „zaplaven“;
  - `obs_conflict`, pokud si odporují (`obs_rp_low > obs_rp_high`).

### 11.6 Spojení s modelovým intervalem

| Případ | Podmínka | Výsledek |
|---|---|---|
| ruční zásah | kompartment má override | pozorování se spočítají, RP se nemění |
| `obs_conflict` | pozorování si odporují | jen příznak, RP se nemění |
| sloučení | `max(rp_low, obs_rp_low) ≤ min(rp_high, obs_rp_high)` | RP se omezí do průniku a přichytí ke třídě uvnitř něj (nejbližší v logaritmu) |
| `model_obs_conflict` | průnik je prázdný: pozorování vyvrací model | použije se jen pozorovaný interval; RP = třída uvnitř něj nejbližší modelové hodnotě |

- **Otevřená strana.** Strana bez pozorování (např. jen testy „vydržel“) zůstane při
  výběru třídy otevřená a pak se nahradí sousední třídou výsledku. To odpovídá
  modelovému intervalu ±1 třída.
  - **Proč.** První verze dosazovala za otevřenou stranu modelovou mez. Při testu
    „vydržel“ nad modelovým intervalem pak mohla vyjít RP *nižší* než pozorování. Chybu
    našel test a je opravená.
- **Žádná třída v intervalu.** Když v úzkém intervalu mezi dvěma třídami žádná třída
  neleží, vezme se nejbližší třída a interval se rozšíří tak, aby ji obsahoval. Vždy tedy
  platí `rp_low ≤ rp ≤ rp_high`; pozorované hodnoty zůstávají v `obs_rp_low` a
  `obs_rp_high`.
- **Opakovatelnost.** Modelové hodnoty se při prvním běhu uloží do sloupců `*_prior` a
  každý další běh vychází z nich. Opakovaný běh s novými pozorováními tak dá totéž co
  první běh (testováno).

**Příklad z testu.** Model dává okruhu 1 hodnotu RP 100 s intervalem [50, 100].
- **Scéna A, „vydržel“ při RP 70:** interval se zúží na [70, 100], RP zůstane 100.
- **Scéna E, „vydržel“ při RP 126:** model je vyvrácen, RP je 200 s intervalem [126, 500].

---

## 12. Validace

`validate.py` srovnává výsledek s referenčními chráněnými územími.

1. **Plochy.**
   - Úplnost (recall) = |výsledek ∩ reference| / |reference| a přesnost (precision) =
     |výsledek ∩ reference| / |výsledek|, obojí na sjednoceních.
   - Pro každý referenční polygon pokrytá plocha a podíl; počty zcela minutých polygonů a
     polygonů pokrytých nad 50 %.
   - Volitelně jen v oblasti `validation_extent`, kde je reference úplná. Mimo ni by
     přesnost neměla smysl.
2. **Doby opakování.**
   - **Párování.** Každý kompartment se srovná s referenčním polygonem, se kterým se
     překrývá nejvíc. Referenční RP se přichytí ke třídám.
   - **Matice záměn** referenční × přiřazená třída.
   - **Podíl v intervalu.** Kompartmenty, jejichž referenční RP leží v `rp_low`…`rp_high`.
   - **Podíl stejné třídy.**
   - **Proč podíl v intervalu.** Hodnotí to, co model skutečně tvrdí: ne bodovou hodnotu,
     ale rozsah.

---

## 13. Výkon a škálování

Všechna čísla pocházejí z tohoto cloudového prostředí: 4 vCPU, 15 GB RAM, syntetická
data. Na vašich strojích se budou lišit hlavně čtení a dekomprese.

| Krok | Naměřeno |
|---|---|
| čtení terénu jednotky 30 M buněk | 8–9 s pásově; 1,7–2,3 s po konverzi A |
| HAND (zahloubení, prolomení, D8) | 13,9 s (dříve 29,3 s) |
| jedno zaplavení (priority flood) | 10,9–11,1 s (dříve 14,6 s) |
| `compartments.py` na jednotku | 56–61 s, špička 2,6–2,7 GB (93–95 B/buňku) |
| `rp_index.py` na jednotku | 1,6–4,6 s, 0,6–0,9 GB |
| `gfm_evidence.py` na jednotku, bez upozornění | asi 5 s + 0,5 s na scénu, 1,3–1,45 GB |
| `gfm_evidence.py` s upozorněními | asi 18 s + 1,6 s na scénu, 2,9 GB |
| `rp_bounds.py`, `validate.py` | sekundy i pro 10 000 kompartmentů |

**Extrapolace.**
- **`compartments.py`.** Čas pracovníka je zhruba 2,4 µs na buňku domény. Při 5 m to dává
  asi 1,6 pracovní minuty na 1000 km² domény. Skript zapisuje extrapolaci z vašich
  testovacích jednotek do `performance.json`.
- **GFM.** Pro 300 jednotek s 20 scénami trvá test kompartmentů asi 75 pracovních minut,
  tedy asi 20 minut na čtyřech pracovnících.
- **Stažení GFM** závisí na službě, kterou jsem odtud nemohl kontaktovat; proto není
  změřené.

---

## 14. Co dokazují testy

Všechny sady prošly při běhu od nuly.

| Sada | Kontroly | Co dokazuje |
|---|---|---|
| `test_kernels` | 16 | uzavřený okruh, mezera, diagonální hráz, vnořené okruhy v obou pořadích výšek, šíření výchozích buněk a pořadí, prolomení a nodata, chyba B1 reprodukovaná bez zahloubení a opravená s ním |
| `test_pipeline` | 40 | celý řetězec fáze A na syntetických zemích: šev mezi jádry (okruh právě jednou), přemostěná a široká mezera, schodovitá hráz, kompartment delší než okno, pravidlo rozporu, opravy B2–B4, řeka celá jako nodata, prázdné hráze |
| `test_stage_b` | 20 | OSM: stavidlo na křížení a čerpací stanice na dolním konci určeném z DTM, takže okruhy jsou chráněné až po kroku OSM; GHSL včetně hledání nejnovější epochy; součty převzorkování v toleranci 1 %; strop, interval, ruční zásah; kvantilové mapování |
| `test_hand_yz` | 7 | výběr dlaždic WorldCover z mřížky; hodnota z land cover shodná s nezávislým výpočtem; přednost Urban Atlasu z dvouvrstvého souboru; otisk vstupů (beze změny se nic nepřepočítá, změna hodnot se pozná); `value_exposed = value_ghsl + value_lulc` |
| `test_stage_c` | 26 | mock STAC (token, stránkování GET i POST); vyřazení skoro suché scény; nodata GFM z konfigurace; výsledky held/flooded/no_test; datumový filtr; upozornění u neznámé suché kapsy a nutnost `hint_exclusion_m`; GloFAS z NetCDF i GRIB se shodným výsledkem; Gumbel; přiřazení podle poměru ploch; všechny případy spojení intervalů; opakovatelnost; validace proti nezávislému výpočtu |

---

## 15. Omezení a otevřené otázky

### 15.1 Omezení metody

- **Smíšené reference HAND.** Buňky stékající k různým řekám se srovnávají v jednom
  zaplavení. U soutoků může drenážní oblast přítoku vypadat chráněná proti hlavní řece;
  to je vlastnost HAND.
- **Useknuté kompartmenty**, které žádné sousední okno nevyvrátí, mohou být uměle
  uzavřené okrajem okna. Jsou označené příznakem a jejich podíl se hlásí.
- **Koruny z DTM.** Protipovodňové zdi užší než rozlišení DTM jsou podhodnocené.
- **Windows.** Fond procesů, měření paměti a extrakce OSM (potřebuje SpatiaLite v GDAL)
  jsou napsané, ale testované jen na Linuxu.
- **Celé říční prvky.** OSM uzavírá celé prvky (8.2).
- **Rozlišení GHSL.** Při 100 m je expozice v malých kompartmentech přibližná.
- **Co radar nevidí.** Voda v hustém lese a zástavbě je vyloučená maskou, takže města
  jsou testována méně často. Voda uvnitř může být i z deště nebo průsaku, takže
  „zaplaven“ může přeceňovat selhání.
- **Rozlišení GloFAS.** Malé řeky (povodí pod několik set km²) GloFAS nerozliší. Jejich
  kompartmenty zůstanou nespárované a ponechají si modelovou RP.
- **Průtok není hladina.** O přelití rozhoduje výška hladiny. Vzdutí od soutoku, ledové
  zátarasy, změny koryta a hystereze vztah mezi RP průtoku a RP hladiny narušují.
- **Čas snímku.** Test „vydržel“ používá denní průměr. Přišla-li kulminace až po přeletu,
  je dolní mez příliš vysoká. Tuto konvenci určuje zadání.

### 15.2 Rozhodnutí k potvrzení

1. **URL zdrojů.** Šablony Geofabrik a GHSL nejsou ověřené. Hostitele byly z vývojového
   prostředí nedostupné.
2. **Služba GFM.** Endpoint STAC, kolekce a názvy vrstev nejsou ověřené. Vyhledávání
   ukázalo `https://stac.eodc.eu/api/v1`, kolekci `GFM` a vrstvy `ensemble_flood_extent`,
   `reference_water_mask`, `ensemble_likelihood`; název masky vyloučení jsem nenašel.
3. **Kódy pixelů GFM** (záplava 1, vyloučení 1, trvalá voda 1, nodata 255) jsou zástupné.
4. **GloFAS.**
   - požadavek z formuláře EWDS;
   - přiřazení dne (`glofas_valid_time_ends_day`);
   - soubor plochy povodí;
   - schválení `cdsapi`.
5. **Pravidla spojení intervalů** (11.6), jak jsem zadání vyložil.
6. **Nový klíč `hint_exclusion_m`** (výchozí 40 m).
7. **Zástupné prahy** `f_*`, rozsah poměru ploch, rozdělení tříd a strop.
8. **Povodeň pro kontrolu v QGIS.** Prosím o určení události.
9. **HAND_YZ.** Převodník tříd na kategorie, hodnoty a křivky kategorií a ověření pole a
   kódů Urban Atlasu; přepínač `lulc.enabled` je ve výchozím stavu vypnutý.

---

## 16. Slovníček

| Česky | Anglicky / v kódu |
|---|---|
| kompartment (chráněné území za hrázemi) | compartment |
| hráz, koruna hráze | levee, crest |
| doba opakování (N-letost) | return period, `rp` |
| úroveň ochrany | `protection_level_m` |
| úroveň přelití | spill level, `spill_undef`, `spill_def` |
| rezerva ochrany | `margin` |
| výška nad nejbližší drenáží | HAND, `hand_undef` |
| výchozí říční buňka | seed |
| vstupní buňka, zdroj | `entry_x/y`, `source_x/y` |
| koridor, jednotka, jádro, lem, okno, doména | corridor, unit, core, halo, window, domain |
| přemostění mezery | bridging, `gap_tolerance_m` |
| stavidlo, vrata, jez, čerpací stanice | sluice gate, lock gate, weir, pumping station |
| zpětné vzdutí | backwater |
| plocha povodí | upstream area, `upstream_area_km2` |
| úplnost, přesnost | recall, precision |
| maska vyloučení, trvalá (referenční) voda | exclusion mask, reference water |
| „vydržel“, „zaplaven“, „neprůkazné“, „bez testu“ | `held`, `flooded`, `inconclusive`, `no_test` |
| upozornění na chybějící bariéru | missing-barrier hint |
| dolní a horní mez | lower and upper bound, `obs_rp_low`, `obs_rp_high` |
| využití a kryt území | land use / land cover (LULC), `value_lulc` |
