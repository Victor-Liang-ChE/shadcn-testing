// Cut these lines from UfcChampionsDisplay.tsx
export interface Champion {
  name: string;
  reignStart: string;
  reignEnd: string | null;
  notes?: string;
}
export interface WeightClassData {
  displayName: string;
  champions: Champion[];
}

export const ufcChampionsData: Record<string, WeightClassData> = {
  /* --------------------------- MEN -------------------------------- */
  heavyweight: {
    displayName: "Heavyweight (206‑265 lbs)",
    champions: [
      { name: "Mark Coleman",   reignStart: "1997-02-07", reignEnd: "1997-07-27" },
      { name: "Maurice Smith",  reignStart: "1997-07-27", reignEnd: "1997-12-21" },
      { name: "Randy Couture",  reignStart: "1997-12-21", reignEnd: "1998-01-15", notes: "stripped (contract)" },
      { name: "Bas Rutten",     reignStart: "1999-05-07", reignEnd: "1999-06-09", notes: "vacated (drop to LHW)" },
      { name: "Kevin Randleman",reignStart: "1999-11-19", reignEnd: "2000-11-17" },
      { name: "Randy Couture",  reignStart: "2000-11-17", reignEnd: "2002-03-22" },
      { name: "Josh Barnett",   reignStart: "2002-03-22", reignEnd: "2002-07-26", notes: "stripped (drug test)" },
      { name: "Ricco Rodriguez",reignStart: "2002-09-27", reignEnd: "2003-02-28" },
      { name: "Tim Sylvia",     reignStart: "2003-02-28", reignEnd: "2003-10-15", notes: "stripped (drug test)" },
      { name: "Frank Mir",      reignStart: "2004-06-19", reignEnd: "2005-08-12", notes: "stripped (injury)" },
      { name: "Andrei Arlovski",reignStart: "2005-08-12", reignEnd: "2006-04-15" },
      { name: "Tim Sylvia",     reignStart: "2006-04-15", reignEnd: "2007-03-03" },
      { name: "Randy Couture",  reignStart: "2007-03-03", reignEnd: "2008-11-15" },
      { name: "Brock Lesnar",   reignStart: "2008-11-15", reignEnd: "2010-10-23" },
      { name: "Cain Velasquez", reignStart: "2010-10-23", reignEnd: "2011-11-12" },
      { name: "Junior dos Santos",reignStart: "2011-11-12",reignEnd: "2012-12-29" },
      { name: "Cain Velasquez", reignStart: "2012-12-29", reignEnd: "2015-06-13" },
      { name: "Fabrício Werdum",reignStart: "2015-06-13", reignEnd: "2016-05-14" },
      { name: "Stipe Miocic",   reignStart: "2016-05-14", reignEnd: "2018-07-07" },
      { name: "Daniel Cormier", reignStart: "2018-07-07", reignEnd: "2019-08-17" },
      { name: "Stipe Miocic",   reignStart: "2019-08-17", reignEnd: "2021-03-27" },
      { name: "Francis Ngannou",reignStart: "2021-03-27", reignEnd: "2023-01-14", notes: "vacated (free agent)" },
      { name: "Jon Jones",      reignStart: "2023-03-04", reignEnd: "Present" },
    ],
  },

  lightheavyweight: {
    displayName: "Light Heavyweight (186‑205 lbs)",
    champions: [
      { name: "Frank Shamrock",   reignStart: "1997-12-21", reignEnd: "1999-11-24" },
      { name: "Tito Ortiz",       reignStart: "2000-04-14", reignEnd: "2003-09-26" },
      { name: "Randy Couture",    reignStart: "2003-09-26", reignEnd: "2004-01-31" },
      { name: "Vitor Belfort",    reignStart: "2004-01-31", reignEnd: "2004-08-21" },
      { name: "Randy Couture",    reignStart: "2004-08-21", reignEnd: "2005-04-16" },
      { name: "Chuck Liddell",    reignStart: "2005-04-16", reignEnd: "2007-05-26" },
      { name: "Quinton Jackson",  reignStart: "2007-05-26", reignEnd: "2008-07-05" },
      { name: "Forrest Griffin",  reignStart: "2008-07-05", reignEnd: "2008-12-27" },
      { name: "Rashad Evans",     reignStart: "2008-12-27", reignEnd: "2009-05-23" },
      { name: "Lyoto Machida",    reignStart: "2009-05-23", reignEnd: "2010-05-08" },
      { name: "Maurício Rua",     reignStart: "2010-05-08", reignEnd: "2011-03-19" },
      { name: "Jon Jones",        reignStart: "2011-03-19", reignEnd: "2015-04-28" },
      { name: "Daniel Cormier",   reignStart: "2015-05-23", reignEnd: "2018-12-28" },
      { name: "Jon Jones",        reignStart: "2018-12-29", reignEnd: "2020-08-17" },
      { name: "Jan Błachowicz",   reignStart: "2020-09-27", reignEnd: "2021-10-30" },
      { name: "Glover Teixeira",  reignStart: "2021-10-30", reignEnd: "2022-06-12" },
      { name: "Jiří Procházka",   reignStart: "2022-06-12", reignEnd: "2022-11-23", notes: "vacated (injury)" },
      { name: "Jamahal Hill",     reignStart: "2023-01-21", reignEnd: "2023-07-14", notes: "vacated (injury)" },
      { name: "Alex Pereira",     reignStart: "2023-11-11", reignEnd: "2025-03-08" },
      { name: "Magomed Ankalaev", reignStart: "2025-03-08", reignEnd: "Present" },
    ],
  },

  middleweight: {
    displayName: "Middleweight (171‑185 lbs)",
    champions: [
      { name: "Dave Menne",        reignStart: "2001-09-28", reignEnd: "2002-01-11" },
      { name: "Murilo Bustamante", reignStart: "2002-01-11", reignEnd: "2002-10-05" },
      { name: "Evan Tanner",       reignStart: "2005-02-05", reignEnd: "2005-06-04" },
      { name: "Rich Franklin",     reignStart: "2005-06-04", reignEnd: "2006-10-14" },
      { name: "Anderson Silva",    reignStart: "2006-10-14", reignEnd: "2013-07-06" },
      { name: "Chris Weidman",     reignStart: "2013-07-06", reignEnd: "2015-12-12" },
      { name: "Luke Rockhold",     reignStart: "2015-12-12", reignEnd: "2016-06-04" },
      { name: "Michael Bisping",   reignStart: "2016-06-04", reignEnd: "2017-11-04" },
      { name: "Georges St-Pierre", reignStart: "2017-11-04", reignEnd: "2017-12-07", notes: "vacated (colitis)" },
      { name: "Robert Whittaker",  reignStart: "2017-12-07", reignEnd: "2019-10-06" },
      { name: "Israel Adesanya",   reignStart: "2019-10-06", reignEnd: "2022-11-12" },
      { name: "Alex Pereira",      reignStart: "2022-11-12", reignEnd: "2023-04-08" },
      { name: "Israel Adesanya",   reignStart: "2023-04-08", reignEnd: "2023-09-10" },
      { name: "Sean Strickland",   reignStart: "2023-09-10", reignEnd: "2024-01-20" },
      { name: "Dricus du Plessis", reignStart: "2024-01-20", reignEnd: "Present" },
    ],
  },

  welterweight: {
    displayName: "Welterweight (156‑170 lbs)",
    champions: [
      { name: "Pat Miletich",  reignStart: "1998-10-16", reignEnd: "2001-05-04" },
      { name: "Carlos Newton", reignStart: "2001-05-04", reignEnd: "2001-11-02" },
      { name: "Matt Hughes",   reignStart: "2001-11-02", reignEnd: "2004-01-31" },
      { name: "B.J. Penn",     reignStart: "2004-01-31", reignEnd: "2004-05-17", notes: "stripped (left UFC)" },
      { name: "Matt Hughes",   reignStart: "2004-10-22", reignEnd: "2006-11-18" },
      { name: "Georges St-Pierre",reignStart: "2006-11-18",reignEnd: "2007-04-07" },
      { name: "Matt Serra",    reignStart: "2007-04-07", reignEnd: "2008-04-19" },
      { name: "Georges St-Pierre",reignStart: "2008-04-19",reignEnd: "2013-12-13", notes: "vacated (hiatus)" },
      { name: "Johny Hendricks",reignStart: "2014-03-15", reignEnd: "2014-12-06" },
      { name: "Robbie Lawler", reignStart: "2014-12-06", reignEnd: "2016-07-30" },
      { name: "Tyron Woodley", reignStart: "2016-07-30", reignEnd: "2019-03-02" },
      { name: "Kamaru Usman",  reignStart: "2019-03-02", reignEnd: "2022-08-20" },
      { name: "Leon Edwards",  reignStart: "2022-08-20", reignEnd: "2024-07-27" },
      { name: "Belal Muhammad",reignStart: "2024-07-27", reignEnd: "2025-05-10" },
      { name: "Jack Della Maddalena",reignStart: "2025-05-10", reignEnd: "Present" },
    ],
  },

  lightweight: {
    displayName: "Lightweight (146‑155 lbs)",
    champions: [
      { name: "Jens Pulver",        reignStart: "2001-02-23", reignEnd: "2002-03-23", notes: "vacated (contract)" },
      { name: "Sean Sherk",         reignStart: "2006-10-14", reignEnd: "2007-12-08", notes: "stripped" },
      { name: "B.J. Penn",          reignStart: "2008-01-19", reignEnd: "2010-04-10" },
      { name: "Frankie Edgar",      reignStart: "2010-04-10", reignEnd: "2012-02-26" },
      { name: "Benson Henderson",   reignStart: "2012-02-26", reignEnd: "2013-08-31" },
      { name: "Anthony Pettis",     reignStart: "2013-08-31", reignEnd: "2015-03-14" },
      { name: "Rafael dos Anjos",   reignStart: "2015-03-14", reignEnd: "2016-07-07" },
      { name: "Eddie Alvarez",      reignStart: "2016-07-07", reignEnd: "2016-11-12" },
      { name: "Conor McGregor",     reignStart: "2016-11-12", reignEnd: "2018-04-07", notes: "stripped (inactivity)" },
      { name: "Khabib Nurmagomedov",reignStart: "2018-04-07", reignEnd: "2021-03-19", notes: "retired" },
      { name: "Charles Oliveira",   reignStart: "2021-05-15", reignEnd: "2022-05-07", notes: "stripped (weight)" },
      { name: "Islam Makhachev",    reignStart: "2022-10-22", reignEnd: "Present" },
    ],
  },

  featherweight: {
    displayName: "Featherweight (136‑145 lbs)",
    champions: [
      { name: "José Aldo",          reignStart: "2010-11-20", reignEnd: "2015-12-12" },
      { name: "Conor McGregor",     reignStart: "2015-12-12", reignEnd: "2016-11-26", notes: "stripped" },
      { name: "José Aldo",          reignStart: "2016-11-26", reignEnd: "2017-06-03" },
      { name: "Max Holloway",       reignStart: "2017-06-03", reignEnd: "2019-12-14" },
      { name: "Alexander Volkanovski",reignStart: "2019-12-14",reignEnd: "2024-02-17" },
      { name: "Ilia Topuria",       reignStart: "2024-02-17", reignEnd: "2025-02-19", notes: "vacated (move to LW)" },
      { name: "Alexander Volkanovski",reignStart: "2025-04-12",reignEnd: "Present" },
    ],
  },

  bantamweight: {
    displayName: "Bantamweight (126‑135 lbs)",
    champions: [
      { name: "Dominick Cruz",     reignStart: "2010-12-16", reignEnd: "2014-01-06" },
      { name: "Renan Barão",       reignStart: "2014-01-06", reignEnd: "2014-05-24" },
      { name: "T.J. Dillashaw",    reignStart: "2014-05-24", reignEnd: "2016-01-17" },
      { name: "Dominick Cruz",     reignStart: "2016-01-17", reignEnd: "2016-12-30" },
      { name: "Cody Garbrandt",    reignStart: "2016-12-30", reignEnd: "2017-11-04" },
      { name: "T.J. Dillashaw",    reignStart: "2017-11-04", reignEnd: "2019-03-20", notes: "vacated (drug test)" },
      { name: "Henry Cejudo",      reignStart: "2019-06-08", reignEnd: "2020-05-24", notes: "vacated (retired)" },
      { name: "Petr Yan",          reignStart: "2020-07-12", reignEnd: "2021-03-06" },
      { name: "Aljamain Sterling", reignStart: "2021-03-06", reignEnd: "2023-08-19" },
      { name: "Sean O'Malley",     reignStart: "2023-08-19", reignEnd: "2024-09-14" },
      { name: "Merab Dvalishvili", reignStart: "2024-09-14", reignEnd: "Present" },
    ],
  },

  flyweight: {
    displayName: "Flyweight (116‑125 lbs)",
    champions: [
      { name: "Demetrious Johnson",reignStart: "2012-09-22", reignEnd: "2018-08-04" },
      { name: "Henry Cejudo",      reignStart: "2018-08-04", reignEnd: "2020-02-29", notes: "vacated" },
      { name: "Deiveson Figueiredo",reignStart: "2020-07-18",reignEnd: "2021-06-12" },
      { name: "Brandon Moreno",    reignStart: "2021-06-12", reignEnd: "2022-01-22" },
      { name: "Deiveson Figueiredo",reignStart: "2022-01-22",reignEnd: "2023-01-21" },
      { name: "Brandon Moreno",    reignStart: "2023-01-21", reignEnd: "2023-07-08" },
      { name: "Alexandre Pantoja", reignStart: "2023-07-08", reignEnd: "Present" },
    ],
  },

  /* -------------------------- WOMEN ------------------------------- */
  women_strawweight: {
    displayName: "Women – Strawweight (106‑115 lbs)",
    champions: [
      { name: "Carla Esparza",   reignStart: "2014-12-12", reignEnd: "2015-03-14" },
      { name: "Joanna Jędrzejczyk",reignStart: "2015-03-14",reignEnd: "2017-11-04" },
      { name: "Rose Namajunas",  reignStart: "2017-11-04", reignEnd: "2019-05-11" },
      { name: "Jéssica Andrade", reignStart: "2019-05-11", reignEnd: "2019-08-31" },
      { name: "Zhang Weili",     reignStart: "2019-08-31", reignEnd: "2021-04-24" },
      { name: "Rose Namajunas",  reignStart: "2021-04-24", reignEnd: "2022-05-07" },
      { name: "Carla Esparza",   reignStart: "2022-05-07", reignEnd: "2022-11-12" },
      { name: "Zhang Weili",     reignStart: "2022-11-12", reignEnd: "Present" },
    ],
  },

  women_flyweight: {
    displayName: "Women – Flyweight (116‑125 lbs)",
    champions: [
      { name: "Nicco Montaño",      reignStart: "2017-12-01", reignEnd: "2018-09-08", notes: "stripped" },
      { name: "Valentina Shevchenko",reignStart: "2018-12-08", reignEnd: "2023-03-04" },
      { name: "Alexa Grasso",       reignStart: "2023-03-04", reignEnd: "2024-09-14" },
      { name: "Valentina Shevchenko",reignStart: "2024-09-14", reignEnd: "Present" },
    ],
  },

  women_bantamweight: {
    displayName: "Women – Bantamweight (126‑135 lbs)",
    champions: [
      { name: "Ronda Rousey",     reignStart: "2012-12-06", reignEnd: "2015-11-15" },
      { name: "Holly Holm",       reignStart: "2015-11-15", reignEnd: "2016-03-05" },
      { name: "Miesha Tate",      reignStart: "2016-03-05", reignEnd: "2016-07-09" },
      { name: "Amanda Nunes",     reignStart: "2016-07-09", reignEnd: "2021-12-11" },
      { name: "Julianna Peña",    reignStart: "2021-12-11", reignEnd: "2022-07-30" },
      { name: "Amanda Nunes",     reignStart: "2022-07-30", reignEnd: "2023-06-20", notes: "retired" },
      /* vacant 2023‑06‑20 → 2024‑01‑20 */
      { name: "Raquel Pennington",reignStart: "2024-01-20", reignEnd: "2024-10-05" },
      { name: "Julianna Peña",    reignStart: "2024-10-05", reignEnd: "Present" },
    ],
  },
};