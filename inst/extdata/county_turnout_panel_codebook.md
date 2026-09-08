# county_turnout_panel — Codebook

`county_turnout_panel` is a tidy panel of U.S. county-level presidential-election
voter turnout and sociodemographic context for the five presidential elections
of **2008, 2012, 2016, 2020, and 2024**. There is one row per county
(county-equivalent) per election year: 15,543 rows, 22 columns.

It generalises the single-election (2016) cross-section `county_turnout` to a
five-wave panel and re-derives every covariate **directly from the U.S. Census
Bureau API** (see *Sources and licensing* below). `turnout` uses the same
votes / CVAP definition as `county_turnout` and reproduces its 2016 values
exactly for 3,103 of 3,107 shared counties (the exceptions are MEDSL vote-return
corrections and this panel's Kansas City handling — see *Known data-quality
flags*).

**Coverage / exclusions**

* 50 states + the District of Columbia. **Alaska is excluded**: MEDSL reports
  Alaska by state-house district, not by borough/census area, so no
  borough-level presidential returns exist. (The original `county_turnout` also
  excludes Alaska.)
* The island territories (American Samoa, Guam, Northern Mariana Islands, Puerto
  Rico, U.S. Virgin Islands) are excluded — no presidential electors / not in
  MEDSL.
* Every covariate is matched to the ACS 5-year vintage whose window **ends
  nearest the election**:

  | election | ACS 5-year | CVAP special tab | USDA ERS RUCC | decennial (rural) |
  |---------:|:-----------|:-----------------|:--------------|:------------------|
  | 2008 | 2006–2010 | 2006–2010 | 2003 | 2010 |
  | 2012 | 2008–2012 | 2008–2012 | 2013 | 2010 |
  | 2016 | 2012–2016 | 2012–2016 | 2013 | 2010 |
  | 2020 | 2016–2020 | 2016–2020 | 2023 | 2020 |
  | 2024 | 2020–2024 | 2020–2024 | 2023 | 2020 |

To regenerate the dataset from scratch, run `Rscript data-raw/county_turnout_panel.R`
from the package root (needs `CENSUS_API_KEY`). See `data-raw/README.md`.

---

## Variables

Listed alphabetically. Each entry gives a description, the year(s)/vintage the
values come from, and the upstream source.

---

**age29andunder_pct**

Description: Percent of the total population aged 29 or younger. Computed from
ACS table **B01001** (Sex by Age) as the sum of the male and female cells for
ages "Under 5" through "25 to 29 years", divided by `B01001_001` × 100.

Year(s): ACS 5-year vintage matched to the election year (see table above).

Source: [U.S. Census Bureau, American Community Survey 5-year estimates, detailed table B01001](https://api.census.gov/data.html) (public domain).

---

**age65andolder_pct**

Description: Percent of the total population aged 65 or older. Computed from ACS
table **B01001** as the sum of the male and female cells for ages "65 and 66
years" through "85 years and over", divided by `B01001_001` × 100.

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed table B01001](https://api.census.gov/data.html) (public domain).

---

**clf_unemploy_pct**

Description: Unemployed persons as a percent of the civilian labor force.
For the 2012–2024 waves: ACS table **B23025** (Employment Status), `B23025_005 /
B23025_003` × 100. For the 2008 wave (2006–2010 ACS), `B23025` does not exist;
the value is built from ACS table **B23001** (Sex by Age by Employment Status) as
Σ(all "…Unemployed" cells) / Σ(all "…Employed" + "…Unemployed" cells) × 100.

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed tables B23025 / B23001](https://api.census.gov/data.html) (public domain).

---

**county**

Description: County (or county-equivalent) name, with a trailing " County" or
" Parish" removed. Independent cities keep their " city" suffix; Alaska boroughs
are not present; Connecticut planning regions (2024 wave only) appear under their
full names (e.g. "Capitol Planning Region"). Not unique on its own — pair with
`fips`.

Year(s): geography of the matched ACS 5-year vintage.

Source: [U.S. Census Bureau, ACS 5-year `NAME` field](https://api.census.gov/data.html) (public domain).

---

**division**

Description: U.S. Census Bureau **division** (9 levels: East North Central, East
South Central, Middle Atlantic, Mountain, New England, Pacific, South Atlantic,
West North Central, West South Central). A fixed state-level classification,
carried over verbatim from the original `county_turnout` dataset.

Year(s): time-invariant.

Source: [U.S. Census Bureau regions and divisions](https://www2.census.gov/geo/pdfs/maps-data/maps/reference/us_regdiv.pdf); inherited from `county_turnout`.

---

**election_year**

Description: Presidential election year — one of 2008, 2012, 2016, 2020, 2024.
The panel key together with `fips`.

Year(s): —

Source: —

---

**female_pct**

Description: Percent of the total population that is female (`B01001_026 /
B01001_001` × 100).

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed table B01001](https://api.census.gov/data.html) (public domain).

---

**fips**

Description: 5-character county FIPS code (state + county), **stored as a
character string with leading zeros** (e.g. `"01001"`). This differs from the
original `county_turnout`, where `fips` is numeric. The panel key together with
`election_year`.

Year(s): geography of the matched ACS 5-year vintage (see *Geographic
harmonisation*).

Source: [U.S. Census Bureau](https://www.census.gov/library/reference/code-lists/ansi.html) (public domain).

---

**foreignborn_pct**

Description: Percent of the total population that is foreign-born (`B05002_013 /
B05002_001` × 100). `B05002_013` ("Foreign born") is a stable line across all
five vintages even though the table's sub-detail changes.

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed table B05002](https://api.census.gov/data.html) (public domain).

---

**geo_group**

Description: A coarse 15-level geographic grouping of states used in the
`pintervals` package examples (e.g. "Eastern Gulf", "Great Lakes Core",
"Pacific Rim"). `pintervals`-specific, not an official classification; carried
over verbatim from the original `county_turnout` dataset.

Year(s): time-invariant.

Source: inherited from `county_turnout` (a `pintervals` addition).

---

**latitude**, **longitude**

Description: County "internal point" (approximate centroid) latitude and
longitude in decimal degrees, rounded to 5 decimals.

Year(s): Census Gazetteer file matched to the election year — 2012 file for the
2008 wave (the Gazetteer series starts in 2012; centroids move < 1 km between
vintages when boundaries are unchanged), 2013 file for 2012, then the same-year
file for 2016 / 2020 / 2024.

Source: [U.S. Census Bureau Gazetteer Files, county `INTPTLAT` / `INTPTLONG`](https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.html) (public domain).

---

**lesscollege_pct**

Description: Percent of the population 25 years and over with **less than a
bachelor's degree** = 100 − Σ(Bachelor's + Master's + Professional + Doctorate) /
total × 100. From ACS table **B15003** (2012–2024 waves) or **B15002**
(2008 wave, where B15003 does not exist — the four degree categories are summed
across the male and female panels).

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed tables B15003 / B15002](https://api.census.gov/data.html) (public domain).

---

**lesshs_pct**

Description: Percent of the population 25 years and over with **less than a high
school diploma** = Σ(all attainment categories below "Regular high school
diploma") / total × 100. From ACS table **B15003** (2012–2024) or **B15002**
(2008 wave; B15002 collapses some grade levels — "Nursery to 4th grade",
"5th and 6th grade", "7th and 8th grade" — which are all counted as below HS).

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed tables B15003 / B15002](https://api.census.gov/data.html) (public domain).

---

**median_hh_inc**

Description: Median household income in the past 12 months, in the inflation-
adjusted dollars of the ACS window's final year (`B19013_001`), stored as an
integer. `NA` where the Census Bureau suppressed the estimate (2 micro-county
rows: Jeff Davis County, TX in 2020; De Baca County, NM in 2024).

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed table B19013](https://api.census.gov/data.html) (public domain).

---

**nonwhite_pct**

Description: Percent of the population that is **not** non-Hispanic White
= 100 − (`B03002_003` / `B03002_001`) × 100, where `B03002_003` is "Not Hispanic
or Latino: White alone".

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed table B03002](https://api.census.gov/data.html) (public domain).

---

**region**

Description: U.S. Census Bureau **region** (4 levels: Midwest, Northeast, South,
West). Fixed state-level classification, carried over verbatim from
`county_turnout`.

Year(s): time-invariant.

Source: [U.S. Census Bureau regions and divisions](https://www2.census.gov/geo/pdfs/maps-data/maps/reference/us_regdiv.pdf); inherited from `county_turnout`.

---

**rural_pct**

Description: Percent of the county population classified as **rural** by the
decennial census "Urban and Rural" table (P2): `P002005 / P002001` × 100 for
2010; `P2_003N / P2_001N` × 100 for 2020. Because the ACS has no urban/rural
classification, this variable is **held fixed from the nearest decennial
census** and is therefore identical across the 2008/2012/2016 waves (2010
census) and across the 2020/2024 waves (2020 census).
`NA` for the nine Connecticut planning regions in 2024: the 2020 census P2 table
is tabulated on the 8 legacy counties, which do not nest within the 9 planning
regions used by the 2020–2024 ACS.

Year(s): 2010 census (2008/2012/2016 waves); 2020 census (2020/2024 waves).

Source: [U.S. Census Bureau, decennial census table P2 (Urban and Rural)](https://api.census.gov/data.html) (public domain).

---

**ruralurban_cc**

Description: USDA Economic Research Service **Rural-Urban Continuum Code**
(integer 1–9; lower = more metropolitan). ERS revises the continuum roughly once
per decade, so vintages are **reused across adjacent elections** (documented, not
hidden):

* **2003** codes → 2008 election
* **2013** codes → 2012 **and** 2016 elections
* **2023** codes → 2020 **and** 2024 elections

Connecticut, 2020 wave: RUCC 2023 is tabulated on the 9 planning regions and does
not join the 8 legacy counties used by the 2016–2020 ACS, so CT's 2020
`ruralurban_cc` is taken from **RUCC 2013** (legacy counties) instead.

Coding:

| code | meaning |
|-----:|:--------|
| 1 | Metro — county in metro area of 1 million population or more |
| 2 | Metro — county in metro area of 250,000 to 1 million population |
| 3 | Metro — county in metro area of fewer than 250,000 population |
| 4 | Nonmetro — urban population ≥ 20,000, adjacent to a metro area |
| 5 | Nonmetro — urban population ≥ 20,000, not adjacent to a metro area |
| 6 | Nonmetro — urban population 2,500–19,999, adjacent to a metro area |
| 7 | Nonmetro — urban population 2,500–19,999, not adjacent to a metro area |
| 8 | Nonmetro — completely rural or < 2,500 urban population, adjacent to a metro area |
| 9 | Nonmetro — completely rural or < 2,500 urban population, not adjacent to a metro area |

Year(s): 2003 / 2013 / 2023 ERS release, as above.

Source: [USDA Economic Research Service, Rural-Urban Continuum Codes](https://www.ers.usda.gov/data-products/rural-urban-continuum-codes/) (U.S. Government work, public domain).

---

**state**

Description: State name (or "District of Columbia").

Year(s): time-invariant.

Source: [U.S. Census Bureau, ACS 5-year `NAME` field](https://api.census.gov/data.html) (public domain).

---

**total_population**

Description: Total county resident population (ACS table **B01003**,
`B01003_001`), stored as an integer.

Year(s): ACS 5-year vintage matched to the election year.

Source: [U.S. Census Bureau, ACS 5-year, detailed table B01003](https://api.census.gov/data.html) (public domain).

---

**turnout**

Description: Voter turnout rate for the presidential election = **total
presidential ballots cast ÷ Citizen Voting-Age Population (CVAP)**.

* Numerator: MEDSL `totalvotes` (the county-wide presidential total, identical
  across all candidate/mode rows for a county-year).
* Denominator: the "Total" line of the Census Bureau CVAP special tabulation
  (`cvap_est`) for the matched ACS 5-year window.

This is the **same definition** as `county_turnout$turnout` (whose own help page
mis-describes the denominator as total population — corrected in `pintervals`
1.3.0). For the 2016 wave, `turnout` here equals `county_turnout$turnout` to
< 0.001 for 3,103 of the 3,107 shared counties. The four larger differences are:
Gila / Graham / Greenlee counties, AZ and Keya Paha County, NE (MEDSL revised its
2016 Arizona/Nebraska county returns after `county_turnout` was built — the panel
uses the corrected v20 returns), plus the four Kansas City-area Missouri counties
(the panel redistributes the Kansas City vote to them; `county_turnout` does
not).

`NA` (14 rows, **retained**): Kalawao County, HI (all 5 waves — no separate
MEDSL returns) and the 9 Connecticut planning regions in 2024 (MEDSL still
reports Connecticut on the 8 legacy counties, which do not match the 2020–2024
ACS/CVAP geography).

**Values > 1 (25 county-years) are dropped.** They occurred only in
micro-counties (e.g. Loving County, TX; Kenedy County, TX; Harding County, NM)
where recorded ballots exceed the noisy 5-year ACS CVAP estimate. These
county-years are removed from the panel by `data-raw/20_assemble_panel.R`
(6 / 3 / 4 / 6 / 6 rows in 2008 / 2012 / 2016 / 2020 / 2024); every non-`NA`
`turnout` in the shipped data is therefore `<= 1`.

Year(s): election-year returns (MEDSL) ÷ matched CVAP window.

Source: [MEDSL, County Presidential Election Returns 2000–2024, `doi:10.7910/DVN/VOQCHQ`](https://doi.org/10.7910/DVN/VOQCHQ) (CC0 1.0); [U.S. Census Bureau CVAP Special Tabulation](https://www.census.gov/programs-surveys/decennial-census/about/voting-rights/cvap.html) (public domain).

---

## Geographic harmonisation

Geography is kept **native to each wave** (the county list of the matched ACS
5-year vintage is the spine). Known boundary/code changes in the 2008–2024
window are handled as follows; nothing is silently dropped.

| change (effective) | handling |
|:-------------------|:---------|
| **Shannon County, SD (46113) → Oglala Lakota County, SD (46102)** — renamed & recoded 2015 | MEDSL (which still uses 46113, even in 2024), RUCC 2013, and the 2010-census rural split are recoded to **46102** for the 2016/2020/2024 waves. The 2008/2012 waves use 46113 (native). |
| **Bedford (independent) city, VA (51515) → merged into Bedford County, VA (51019)**, 2013 | 2008 & 2012 waves: 51515 is a separate row (native, pre-merge). 2016 wave: MEDSL votes, RUCC, and the census rural counts for 51515 are folded into **51019**. 2020/2024: 51019 only. |
| **Wade Hampton Census Area, AK (02270) → Kusilvak Census Area, AK (02158)**, 2015 | Alaska is excluded from the panel; noted for completeness only. |
| **Valdez-Cordova Census Area, AK (02261) split → Chugach (02063) + Copper River (02066)**, 2019 | Alaska is excluded from the panel; noted for completeness only. |
| **Connecticut: 8 counties → 9 COG-based planning regions**, 2022-vintage products | Kept native per wave: **2008–2020 waves use the 8 legacy counties (09001–09015); the 2024 wave uses the 9 planning regions (09110–09190).** Connecticut is therefore **not** a consistent 5-wave panel. Consequences: 2024 `turnout` is `NA` for the 9 planning regions (MEDSL still reports the 8 counties); 2024 `rural_pct` is `NA` (2020 census on legacy counties); 2020 `ruralurban_cc` comes from RUCC 2013. |
| **Kansas City, MO** reports as its own unit (FIPS `2938000`; mis-coded `36000` in the 2024 MEDSL file — an upstream error) straddling Cass/Clay/Jackson/Platte counties, whose own rows exclude the Kansas City vote | Kansas City's presidential ballots are redistributed to the four counties in proportion to each county's share of the Kansas City *place* population from the nearest decennial census, then added to those counties' numerators. |

Per-year county counts (investigate, do not assume constant). Native geography
per wave, then the `turnout > 1` county-years dropped by
`data-raw/20_assemble_panel.R`:

| year | native geo | − turnout > 1 | shipped rows |
|-----:|-----------:|--------------:|-------------:|
| 2008 | 3114 | 6 | 3108 |
| 2012 | 3114 | 3 | 3111 |
| 2016 | 3113 (−1: Bedford city → Bedford County) | 4 | 3109 |
| 2020 | 3113 | 6 | 3107 |
| 2024 | 3114 (+1: Connecticut 8 counties → 9 planning regions) | 6 | 3108 |

## Known data-quality flags

* **`turnout` > 1** in 25 micro-county-years — **dropped** from the panel; see
  the `turnout` entry.
* **`turnout` is `NA`** for 14 rows — Kalawao Co., HI and CT-2024 planning
  regions; retained. See the `turnout` entry.
* **`median_hh_inc` is `NA`** for 2 micro-county rows (Census-suppressed).
* **`rural_pct` is `NA`** for the 9 CT-2024 planning regions.
* **Small-county ACS noise**: 5-year ACS estimates for counties with a few
  hundred residents (e.g. Loving Co., TX; Kalawao Co., HI) have very large
  margins of error. Extreme values such as `age65andolder_pct ≈ 82` (Loving Co.)
  are genuine point estimates, not errors.
* **National CVAP turnout** implied by the panel (CVAP-weighted mean of
  `turnout`, `turnout > 1` rows already removed): 2008 ≈ 0.62, 2012 ≈ 0.60,
  2016 ≈ 0.61, 2020 ≈ 0.68, 2024 ≈ 0.64 — within ~1–2 points of published
  national figures (the CVAP window is centred a year or two before each
  election, slightly depressing the denominator).

## Sources and licensing

| source | used for | license |
|:-------|:---------|:--------|
| MEDSL, *County Presidential Election Returns 2000–2024*, Harvard Dataverse `doi:10.7910/DVN/VOQCHQ` (v20.0) | `turnout` numerator | **CC0 1.0** (verified via the Dataverse API) |
| U.S. Census Bureau, CVAP Special Tabulation from the ACS (`census.gov/rdo`) | `turnout` denominator | U.S. Government work — **public domain** |
| U.S. Census Bureau, ACS 5-year **detailed tables**, pulled **directly from the Census Bureau API** | all demographic covariates | U.S. Government work — **public domain** |
| U.S. Census Bureau, decennial census table P2 (Census API) | `rural_pct` | U.S. Government work — **public domain** |
| U.S. Census Bureau, Gazetteer Files | `longitude`, `latitude` | U.S. Government work — **public domain** |
| USDA Economic Research Service, Rural-Urban Continuum Codes | `ruralurban_cc` | U.S. Government work — **public domain** |
| original `pintervals::county_turnout` | `division`, `region`, `geo_group` | see that dataset |

**Why the ACS covariates are pulled from the Census Bureau API and not IPUMS
NHGIS:** NHGIS's terms of use prohibit redistributing extracted data outside
narrow journal-specific exceptions, so ACS data obtained through NHGIS cannot be
bundled in a CRAN package. The identical ACS detailed tables obtained directly
from the Census Bureau's own API are a U.S. Government work in the public domain
and carry no such restriction. The reconstruction was validated against the
original `county_turnout` (built from the MEDSL 2018 context file, which used
NHGIS): for the 2012–2016 ACS vintage, all 10 shared covariates match to within
floating-point error across all 3,107 shared counties.
