/******************************************************************************/

size_t CodonSiteTools::numberOfSubsitutions(const Site& site, double freqmin)
throw (AlphabetException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::numberOfSubsitutions: alphabet is not CodonAlphabet", site.getAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::numberOfSubsitutions: Incorrect specified site", &site);

  if (SiteTools::isConstant(site))
    return 0;
  Site* newsite;
  if (freqmin > 1. / static_cast<double>(site.size()))
    newsite = CodonSiteTools::generateCodonSiteWithoutRareVariant(site, freqmin);
  else
    newsite = new Site(site);
  // Computation
  if (SiteTools::hasGap(*newsite))
    return 0;
  vector<int> pos1, pos2, pos3;

  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(site.getAlphabet());

  for (size_t i = 0; i < newsite->size(); i++)
  {
    pos1.push_back(ca->getFirstPosition(newsite->getValue(i)));
    pos2.push_back(ca->getSecondPosition(newsite->getValue(i)));
    pos3.push_back(ca->getThirdPosition(newsite->getValue(i)));
  }

  const NucleicAlphabet* na = ca->getNucleicAlphabet();

  Site s1(pos1, na), s2(pos2, na), s3(pos3, na);
  size_t Scodon = SiteTools::getNumberOfDistinctCharacters(*newsite) - 1;
  size_t Sbase = SiteTools::getNumberOfDistinctCharacters(s1) + SiteTools::getNumberOfDistinctCharacters(s2) + SiteTools::getNumberOfDistinctCharacters(s3) - 3;
  delete newsite;
  if (Scodon >= Sbase)
    return Scodon;
  else
    return Sbase;
}

/******************************************************************************/

size_t CodonSiteTools::numberOfNonSynonymousSubstitutions(const Site& site, const GeneticCode& gc, double freqmin)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::numberOfNonSynonymousSubstitutions: alphabet is not CodonAlphabet", site.getAlphabet());
  if (typeid(site.getAlphabet()) != typeid(gc.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::numberOfNonSynonymousSubstitutions: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gc.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::numberOfNonSynonymousSubstitutions: Incorrect specified site", &site);

  if (SiteTools::isConstant(site))
    return 0;
  Site* newsite;
  if (freqmin > 1. / static_cast<double>(site.size()))
    newsite = generateCodonSiteWithoutRareVariant(site, freqmin);
  else
    newsite = new Site(site);
  if (SiteTools::hasGap(*newsite))
    return 0;
  // computation
  map<int, size_t> count;
  SiteTools::getCounts(*newsite, count);
  size_t NaSup = 0;
  size_t Nminmin = 10;

  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(site.getAlphabet());

  for (map<int, size_t>::iterator it1 = count.begin(); it1 != count.end(); it1++)
  {
    size_t Nmin = 10;
    for (map<int, size_t>::iterator it2 = count.begin(); it2 != count.end(); it2++)
    {
      size_t Ntot = numberOfDifferences(it1->first, it2->first, *ca);
      size_t Ns = (size_t)numberOfSynonymousDifferences(it1->first, it2->first, gc, true);
      if (Nmin > Ntot - Ns && it1->first != it2->first)
        Nmin = Ntot - Ns;
    }
    NaSup += Nmin;
    if (Nmin < Nminmin)
      Nminmin = Nmin;
  }
  delete newsite;
  return NaSup - Nminmin;
}

/******************************************************************************/

size_t SiteTools::getNumberOfDistinctCharacters(const Site& site) throw (EmptySiteException)
{
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("SiteTools::getNumberOfDistinctCharacters(): Incorrect specified site, size must be > 0", &site);
  // For all site's characters
  if (SiteTools::isConstant(site))
    return 1;
  map<int, size_t> counts;
  SymbolListTools::getCounts(site, counts);
  int s = 0;
  for (map<int, size_t>::iterator it = counts.begin(); it != counts.end(); it++)
  {
    if (it->second != 0)
      s++;
  }
  return s;
}


