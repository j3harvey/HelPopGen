/******************************************************************************/

size_t CodonSiteTools::numberOfSubsitutions(const Site& site, double freqmin)
{
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
{
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


// ******************************************************************************

size_t SequenceStatistics::polymorphicSiteNumber(const PolymorphismSequenceContainer& psc, bool gapflag, bool ignoreUnknown)
{
  size_t S = 0;
  const Site* site = 0;
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (!SiteTools::isConstant(*site, ignoreUnknown))
    {
      S++;
    }
  }
  delete si;
  return S;
}

size_t SequenceStatistics::parsimonyInformativeSiteNumber(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  size_t S = 0;
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (SiteTools::isParsimonyInformativeSite(*site))
    {
      S++;
    }
  }
  delete si;
  return S;
}

size_t SequenceStatistics::countSingleton(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  size_t nus = 0;
  const Site* site = 0;
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    nus += getSingletonNumber_(*site);
  }
  delete si;
  return nus;
}

size_t SequenceStatistics::tripletNumber(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  int S = 0;
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (SiteTools::isTriplet(*site))
    {
      S++;
    }
  }

  delete si;
  return S;
}

size_t SequenceStatistics::totNumberMutations(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  size_t tnm = 0;
  const Site* site = 0;
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    tnm += getMutationNumber_(*site);
  }
  delete si;
  return tnm;
}

size_t SequenceStatistics::totMutationsExternalBranchs(
  const PolymorphismSequenceContainer& ing,
  const PolymorphismSequenceContainer& outg) throw (Exception)
{
  if (ing.getNumberOfSites() != outg.getNumberOfSites())
    throw Exception("ing and outg must have the same size");
  size_t nmuts = 0;
  const Site* site_in = 0;
  const Site* site_out = 0;
  ConstSiteIterator* si = 0;
  ConstSiteIterator* so = 0;
  si = new SimpleSiteContainerIterator(ing);
  so = new SimpleSiteContainerIterator(outg);
  while (si->hasMoreSites())
  {
    site_in = si->nextSite();
    site_out = so->nextSite();
    // use fully resolved sites
    if (SiteTools::isComplete(*site_in) &&  SiteTools::isComplete(*site_out))
      nmuts += getDerivedSingletonNumber_(*site_in, *site_out);                                                                   // singletons that are not in outgroup
  }
  delete si;
  delete so;
  return nmuts;
}

double SequenceStatistics::heterozygosity(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  ConstSiteIterator* si = 0;
  const Site* site = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  double S = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    S += SiteTools::heterozygosity(*site);
  }
  delete si;
  return S;
}

double SequenceStatistics::squaredHeterozygosity(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  ConstSiteIterator* si = 0;
  const Site* site = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  double S = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    double h = SiteTools::heterozygosity(*site);
    S += h * h;
  }
  delete si;
  return S;
}

// ******************************************************************************
// Synonymous and non-synonymous polymorphism
// ******************************************************************************

size_t SequenceStatistics::stopCodonSiteNumber(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  /*
   * Sylvain Gaillard 17/03/2010
   * What if the Alphabet is not a codon alphabet?
   */
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new NoGapSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  size_t S = 0;
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (CodonSiteTools::hasStop(*site))
      S++;
  }
  delete si;
  return S;
}

size_t SequenceStatistics::monoSitePolymorphicCodonNumber(const PolymorphismSequenceContainer& psc, bool stopflag, bool gapflag)
{
  ConstSiteIterator* si = 0;
  if (stopflag)
    si = new CompleteSiteContainerIterator(psc);
  else
  {
    if (gapflag)
      si = new NoGapSiteContainerIterator(psc);
    else
      si = new SimpleSiteContainerIterator(psc);
  }
  size_t S = 0;
  const Site* site;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (CodonSiteTools::isMonoSitePolymorphic(*site))
      S++;
  }
  delete si;
  return S;
}

size_t SequenceStatistics::synonymousPolymorphicCodonNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc)
{
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  size_t S = 0;
  const Site* site;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (CodonSiteTools::isSynonymousPolymorphic(*site, gc))
      S++;
  }
  delete si;
  return S;
}

double SequenceStatistics::watterson75Synonymous(const PolymorphismSequenceContainer& psc, const GeneticCode& gc)
{
  double ThetaW = 0.;
  size_t n = psc.getNumberOfSequences();
  size_t S = synonymousSubstitutionsNumber(psc, gc);
  map<string, double> values = getUsefullValues_(n);
  ThetaW = static_cast<double>(S) / values["a1"];
  return ThetaW;
}

double SequenceStatistics::watterson75NonSynonymous(const PolymorphismSequenceContainer& psc, const GeneticCode& gc)
{
  double ThetaW;
  size_t n = psc.getNumberOfSequences();
  size_t S = nonSynonymousSubstitutionsNumber(psc, gc);
  map<string, double> values = getUsefullValues_(n);
  ThetaW = static_cast<double>(S) / values["a1"];
  return ThetaW;
}

double SequenceStatistics::piSynonymous(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, bool minchange)
{
  double S = 0.;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    S += CodonSiteTools::piSynonymous(*site, gc, minchange);
  }
  delete si;
  return S;
}

double SequenceStatistics::piNonSynonymous(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, bool minchange)
{
  double S = 0.;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    S += CodonSiteTools::piNonSynonymous(*site, gc, minchange);
  }
  delete si;
  return S;
}

double SequenceStatistics::meanSynonymousSitesNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double ratio)
{
  double S = 0.;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    S += CodonSiteTools::meanNumberOfSynonymousPositions(*site, gc, ratio);
  }
  delete si;
  return S;
}

double SequenceStatistics::meanNonSynonymousSitesNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double ratio)
{
  double S = 0.;
  int n = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    n = n + 3;
    S += CodonSiteTools::meanNumberOfSynonymousPositions(*site, gc, ratio);
  }
  delete si;
  return static_cast<double>(n - S);
}

size_t SequenceStatistics::synonymousSubstitutionsNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double freqmin)
{
  size_t St = 0, Sns = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    St += CodonSiteTools::numberOfSubsitutions(*site, freqmin);
    Sns += CodonSiteTools::numberOfNonSynonymousSubstitutions(*site, gc, freqmin);
  }
  delete si;
  return St - Sns;
}

size_t SequenceStatistics::nonSynonymousSubstitutionsNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double freqmin)
{
  size_t Sns = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    Sns += CodonSiteTools::numberOfNonSynonymousSubstitutions(*site, gc, freqmin);
  }
  delete si;
  return Sns;
}

vector<size_t> SequenceStatistics::fixedDifferences(const PolymorphismSequenceContainer& pscin, const PolymorphismSequenceContainer& pscout, PolymorphismSequenceContainer& psccons, const GeneticCode& gc)
{
  ConstSiteIterator* siIn = new CompleteSiteContainerIterator(pscin);
  ConstSiteIterator* siOut = new CompleteSiteContainerIterator(pscout);
  ConstSiteIterator* siCons = new CompleteSiteContainerIterator(psccons);
  const Site* siteIn = 0;
  const Site* siteOut = 0;
  const Site* siteCons = 0;
  size_t NfixS = 0;
  size_t NfixA = 0;
  while (siIn->hasMoreSites())
  {
    siteIn = siIn->nextSite();
    siteOut = siOut->nextSite();
    siteCons = siCons->nextSite();
    vector<size_t> v = CodonSiteTools::fixedDifferences(*siteIn, *siteOut, siteCons->getValue(0), siteCons->getValue(1), gc);
    NfixS += v[0];
    NfixA += v[1];
  }
  vector<size_t> v(2);
  v[0] = NfixS;
  v[1] = NfixA;
  delete siIn;
  delete siOut;
  delete siCons;
  return v;
}

vector<size_t> SequenceStatistics::MKtable(const PolymorphismSequenceContainer& ingroup, const PolymorphismSequenceContainer& outgroup, const GeneticCode& gc, double freqmin)
{
  PolymorphismSequenceContainer psctot(ingroup);
  for (size_t i = 0; i < outgroup.getNumberOfSequences(); i++)
  {
    psctot.addSequence(outgroup.getSequence(i));
    psctot.setAsOutgroupMember(i + ingroup.getNumberOfSequences());
  }
  const PolymorphismSequenceContainer* psccomplet = PolymorphismSequenceContainerTools::getCompleteSites(psctot);
  const PolymorphismSequenceContainer* pscin = PolymorphismSequenceContainerTools::extractIngroup(*psccomplet);
  const PolymorphismSequenceContainer* pscout = PolymorphismSequenceContainerTools::extractOutgroup(*psccomplet);
  const Sequence* consensusIn = SiteContainerTools::getConsensus(*pscin, "consensusIn");
  const Sequence* consensusOut = SiteContainerTools::getConsensus(*pscout, "consensusOut");
  PolymorphismSequenceContainer* consensus = new PolymorphismSequenceContainer(ingroup.getAlphabet());
  consensus->addSequence(*consensusIn);
  consensus->addSequence(*consensusOut);
  vector<size_t> u = SequenceStatistics::fixedDifferences(*pscin, *pscout, *consensus, gc);
  vector<size_t> v(4);
  v[0] = SequenceStatistics::nonSynonymousSubstitutionsNumber(*pscin, gc, freqmin);
  v[1] = SequenceStatistics::synonymousSubstitutionsNumber(*pscin, gc, freqmin);
  v[2] = u[1];
  v[3] = u[0];
  delete consensus;
  if (psccomplet)
  {
    delete psccomplet;
  }
  if (pscin)
  {
    delete pscin;
  }
  if (pscout)
  {
    delete pscout;
  }
  if (consensusIn)
  {
    delete consensusIn;
  }
  if (consensusOut)
  {
    delete consensusOut;
  }
  return v;
}

double SequenceStatistics::neutralityIndex(const PolymorphismSequenceContainer& ingroup, const PolymorphismSequenceContainer& outgroup, const GeneticCode& gc, double freqmin)
{
  vector<size_t> v = SequenceStatistics::MKtable(ingroup, outgroup, gc, freqmin);
  if (v[1] != 0 && v[2] != 0)
    return static_cast<double>(v[0] * v[3]) / static_cast<double>(v[1] * v[2]);
  else
    return -1;
}

/******************************************************************************/

size_t CodonSiteTools::numberOfDifferences(int i, int j, const CodonAlphabet& ca)
{
  size_t nbdif = 0;
  if (ca.getFirstPosition(i) != ca.getFirstPosition(j))
    nbdif++;
  if (ca.getSecondPosition(i) != ca.getSecondPosition(j))
    nbdif++;
  if (ca.getThirdPosition(i) != ca.getThirdPosition(j))
    nbdif++;
  return nbdif;
}

/******************************************************************************/

double CodonSiteTools::numberOfSynonymousDifferences(int i, int j, const GeneticCode& gc, bool minchange)
{
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet());
  vector<int> ci = ca->getPositions(i);
  vector<int> cj = ca->getPositions(j);

  switch (numberOfDifferences(i, j, *ca))
  {
  case 0: return 0;
  case 1:
  {
    if (gc.areSynonymous(i, j))
      return 1;
    return 0;
  }
  case 2:
  {
    if (gc.areSynonymous(i, j))
      return 2;
    vector<double> path(2, 0); // Vector of number of synonymous changes per path (2 here)
    vector<double> weight(2, 1); // Weight to exclude path through stop codon
    if (ci[0] == cj[0])
    {
      int trans1 = ca->getCodon(ci[0], cj[1], ci[2]); // transitory codon between NcNiNi et NcNjNj: NcNjNi, Nc = identical site
      int trans2 = ca->getCodon(ci[0], ci[1], cj[2]); // transitory codon between NcNiNi et NcNjNj: NcNiNj, Nc = identical site
      if (!ca->isStop(trans1))
      {
        if (gc.areSynonymous(i, trans1))
          path[0]++;
        if (gc.areSynonymous(trans1, j))
          path[0]++;
      }
      else
        weight[0] = 0;
      if (!ca->isStop(trans2))
      {
        if (gc.areSynonymous(i, trans2))
          path[1]++;
        if (gc.areSynonymous(trans2, j))
          path[1]++;
      }
      else
        weight[1] = 0;
    }
    if (ci[1] == cj[1])
    {
      int trans1 = ca->getCodon(cj[0], ci[1], ci[2]); // transitory codon between NiNcNi et NjNcNj: NjNcNi, Nc = identical site
      int trans2 = ca->getCodon(ci[0], ci[1], cj[2]); // transitory codon between NiNcNi et NjNcNj: NiNcNj, Nc = identical site
      if (!ca->isStop(trans1))
      {
        if (gc.areSynonymous(i, trans1))
          path[0]++;
        if (gc.areSynonymous(trans1, j))
          path[0]++;
      }
      else
        weight[0] = 0;
      if (!ca->isStop(trans2))
      {
        if (gc.areSynonymous(i, trans2))
          path[1]++;
        if (gc.areSynonymous(trans2, j))
          path[1]++;
      }
      else
        weight[1] = 0;
    }
    if (ci[2] == cj[2])
    {
      int trans1 = ca->getCodon(cj[0], ci[1], ci[2]); // transitory codon between NiNiNc et NjNjNc: NjNiNc, Nc = identical site
      int trans2 = ca->getCodon(ci[0], cj[1], ci[2]); // transitory codon between NiNiNc et NjNjNc: NiNjNc, Nc = identical site
      if (!ca->isStop(trans1))
      {
        if (gc.areSynonymous(i, trans1))
          path[0]++;
        if (gc.areSynonymous(trans1, j))
          path[0]++;
      }
      else
        weight[0] = 0;
      if (!ca->isStop(trans2))
      {
        if (gc.areSynonymous(i, trans2))
          path[1]++;
        if (gc.areSynonymous(trans2, j))
          path[1]++;
      }
      else
        weight[1] = 0;
    }
    if (minchange)
      return VectorTools::max(path);

    double nbdif = 0;
    for (size_t k = 0; k < 2; k++)
    {
      nbdif += path[k] * weight[k];
    }

    return nbdif / VectorTools::sum(weight);
  }
  case 3:
  {
    vector<double> path(6, 0);
    vector<double> weight(6, 1);
    // First transitory codons
    int trans100 = ca->getCodon(cj[0], ci[1], ci[2]);
    int trans010 = ca->getCodon(ci[0], cj[1], ci[2]);
    int trans001 = ca->getCodon(ci[0], ci[1], cj[2]);
    // Second transitory codons
    int trans110 = ca->getCodon(cj[0], cj[1], ci[2]);
    int trans101 = ca->getCodon(cj[0], ci[1], cj[2]);
    int trans011 = ca->getCodon(ci[0], cj[1], cj[2]);
    // Paths
    if (!ca->isStop(trans100))
    {
      if (gc.areSynonymous(i, trans100))
      {
        path[0]++; path[1]++;
      }
      if (!ca->isStop(trans110))
      {
        if (gc.areSynonymous(trans100, trans110))
          path[0]++;
        if (gc.areSynonymous(trans110, j))
          path[0]++;
      }
      else
        weight[0] = 0;
      if (!ca->isStop(trans101))
      {
        if (gc.areSynonymous(trans100, trans101))
          path[1]++;
        if (gc.areSynonymous(trans101, j))
          path[1]++;
      }
      else
        weight[1] = 0;
    }
    else
    {
      weight[0] = 0; weight[1] = 0;
    }
    if (!ca->isStop(trans010))
    {
      if (gc.areSynonymous(i, trans010))
      {
        path[2]++; path[3]++;
      }
      if (!ca->isStop(trans110))
      {
        if (gc.areSynonymous(trans010, trans110))
          path[2]++;
        if (gc.areSynonymous(trans110, j))
          path[2]++;
      }
      else
        weight[2] = 0;
      if (!ca->isStop(trans011))
      {
        if (gc.areSynonymous(trans010, trans011))
          path[3]++;
        if (gc.areSynonymous(trans011, j))
          path[3]++;
      }
      else
        weight[3] = 0;
    }
    else
    {
      weight[2] = 0; weight[3] = 0;
    }
    if (!ca->isStop(trans001))
    {
      if (gc.areSynonymous(i, trans001))
      {
        path[4]++; path[5]++;
      }
      if (!ca->isStop(trans101))
      {
        if (gc.areSynonymous(trans001, trans101))
          path[4]++;
        if (gc.areSynonymous(trans101, j))
          path[4]++;
      }
      else
        weight[4] = 0;
      if (!ca->isStop(trans011))
      {
        if (gc.areSynonymous(trans001, trans011))
          path[5]++;
        if (gc.areSynonymous(trans011, j))
          path[5]++;
      }
      else
        weight[5] = 0;
    }
    else
    {
      weight[4] = 0; weight[5] = 0;
    }
    if (minchange)
      return VectorTools::max(path);

    double nbdif = 0;
    for (size_t k = 0; k < 6; k++)
    {
      nbdif += path[k] * weight[k];
    }

    return nbdif / VectorTools::sum(weight);
  }
  }
  // This line is never reached but sends a warning if not there:
  return 0.;
}

/******************************************************************************/

double CodonSiteTools::piSynonymous(const Site& site, const GeneticCode& gc, bool minchange)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::piSynonymous: alphabet is not CodonAlphabet", site.getAlphabet());
  if (typeid(site.getAlphabet()) != typeid(gc.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::piSynonymous: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gc.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::piSynonymous: Incorrect specified site", &site);

  // General polymorphism checking
  if (SiteTools::isConstant(site))
    return 0;
  // Computation
  map<int, double> freq;
  SiteTools::getFrequencies(site, freq);
  double pi = 0;
  for (map<int, double>::iterator it1 = freq.begin(); it1 != freq.end(); it1++)
  {
    for (map<int, double>::iterator it2 = freq.begin(); it2 != freq.end(); it2++)
    {
      pi += (it1->second) * (it2->second) * (numberOfSynonymousDifferences(it1->first, it2->first, gc, minchange));
    }
  }
  size_t n = site.size();
  return pi * static_cast<double>(n / (n - 1));
}

/******************************************************************************/

double CodonSiteTools::piNonSynonymous(const Site& site, const GeneticCode& gc, bool minchange)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::piNonSynonymous: alphabet is not CodonAlphabet", site.getAlphabet());
  if (typeid(site.getAlphabet()) != typeid(gc.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::piNonSynonymous: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gc.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::piSynonymous: Incorrect specified site", &site);

  // General polymorphism checking
  if (SiteTools::isConstant(site))
    return 0;
  if (isSynonymousPolymorphic(site, gc))
    return 0;
  // Computation
  map<int, double> freq;
  SiteTools::getFrequencies(site, freq);
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(site.getAlphabet());
  double pi = 0;
  for (map<int, double>::iterator it1 = freq.begin(); it1 != freq.end(); it1++)
  {
    for (map<int, double>::iterator it2 = freq.begin(); it2 != freq.end(); it2++)
    {
      double nbtot = static_cast<double>(numberOfDifferences(it1->first, it2->first, *ca));
      double nbsyn = numberOfSynonymousDifferences(it1->first, it2->first, gc, minchange);
      pi += (it1->second) * (it2->second) * (nbtot - nbsyn);
    }
  }
  size_t n = site.size();
  return pi * static_cast<double>(n / (n - 1));
}

/******************************************************************************/

double CodonSiteTools::numberOfSynonymousPositions(int i, const GeneticCode& gc, double ratio) throw (Exception)
{
  try
  {
    const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet());
    if (ca->isStop(i))
      return 0;
    if (ca->isUnresolved(i))
      return 0;
    double nbsynpos = 0.0;
    vector<int> codon = ca->getPositions(i);
    int acid = gc.translate(i);
    for (int pos = 0; pos < 3; pos++)
    {
      for (int an = 0; an < 4; an++)
      {
        if (an == codon[pos])
          continue;
        vector<int> mutcodon = codon;
        mutcodon[pos] = an;
        int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);
        if (ca->isStop(intcodon))
          continue;
        int altacid = gc.translate(intcodon);
        if (altacid == acid)   // if synonymous
        {
          if (((codon[pos] == 0 || codon[pos] == 2) && (mutcodon[pos] == 1 || mutcodon[pos] == 3)) ||
              ((codon[pos] == 1 || codon[pos] == 3) && (mutcodon[pos] == 0 || mutcodon[pos] == 2)))   // if it is a transversion
          {
            nbsynpos = nbsynpos + 1 / (ratio + 2);
          }
          else   // if transition
          {
            nbsynpos = nbsynpos + ratio / (ratio + 2);
          }
        }
      }
    }
    return nbsynpos;
  }
  catch (...)
  {} // !!!!! en cas d'exception, plante! il faudrait forwarder l'exception
  // This line is never reached but sends a warning if not there:
  return 0.;
}

/******************************************************************************/

double CodonSiteTools::meanNumberOfSynonymousPositions(const Site& site, const GeneticCode& gc, double ratio)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::meanNumberOfSynonymousPositions: alphabet is not CodonAlphabet", site.getAlphabet());
  if (typeid(site.getAlphabet()) != typeid(gc.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::meanNumberOfSynonymousPositions: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gc.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::meanNumberOfSynonymousPositions: Incorrect specified site", &site);

  // Computation
  double NbSyn = 0;
  map<int, double> freq;
  SiteTools::getFrequencies(site, freq);
  for (map<int, double>::iterator it = freq.begin(); it != freq.end(); it++)
  {
    NbSyn += (it->second) * numberOfSynonymousPositions(it->first, gc, ratio);
  }
  return NbSyn;
}

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

vector<size_t> CodonSiteTools::fixedDifferences(const Site& siteIn, const Site& siteOut, int i, int j, const GeneticCode& gc)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(siteIn.getAlphabet()))
    throw AlphabetException("CodonSiteTools::fixedDifferences: alphabet is not CodonAlphabet (siteIn)", siteIn.getAlphabet());
  if (!AlphabetTools::isCodonAlphabet(siteOut.getAlphabet()))
    throw AlphabetException("CodonSiteTools::fixedDifferences: alphabet is not CodonAlphabet (siteOut)", siteOut.getAlphabet());
  if (typeid(siteIn.getAlphabet()) != typeid(gc.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::fixedDifferences: siteIn and genetic code have not the same codon alphabet.", siteIn.getAlphabet(), gc.getSourceAlphabet());
  if (typeid(siteOut.getAlphabet()) != typeid(gc.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::fixedDifferences: siteOut and genetic code have not the same codon alphabet.", siteOut.getAlphabet(), gc.getSourceAlphabet());
  // Empty site checking
  if (siteIn.size() == 0)
    throw EmptySiteException("CodonSiteTools::getFixedDifferences Incorrect specified site", &siteIn);
  if (siteOut.size() == 0)
    throw EmptySiteException("CodonSiteTools::getFixedDifferences Incorrect specified site", &siteOut);

  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet());

  size_t Ntot = numberOfDifferences(i, j, *ca);
  size_t Ns = (size_t) numberOfSynonymousDifferences(i, j, gc, true);
  size_t Na = Ntot - Ns;
  size_t Nfix = Ntot;
  vector<int> pos1in, pos2in, pos3in, pos1out, pos2out, pos3out;

  for (size_t k = 0; k < siteIn.size(); k++)
  {
    pos1in.push_back(ca->getFirstPosition(siteIn[k]));
    pos2in.push_back(ca->getSecondPosition(siteIn[k]));
    pos3in.push_back(ca->getThirdPosition(siteIn[k]));
    pos1out.push_back(ca->getFirstPosition(siteOut[k]));
    pos2out.push_back(ca->getSecondPosition(siteOut[k]));
    pos3out.push_back(ca->getThirdPosition(siteOut[k]));
  }
  const NucleicAlphabet* na = ca->getNucleicAlphabet();

  Site s1in(pos1in, na), s2in(pos2in, na), s3in(pos3in, na);
  Site s1out(pos1out, na), s2out(pos2out, na), s3out(pos3out, na);
  bool test1 = false;
  bool test2 = false;
  bool test3 = false;
  if ( (!SiteTools::isConstant(s1in) || !SiteTools::isConstant(s1out)) && ca->getFirstPosition(i) != ca->getFirstPosition(j) )
  {
    test1 = true;
    Nfix--;
  }
  if ( (!SiteTools::isConstant(s2in) || !SiteTools::isConstant(s2out)) && ca->getSecondPosition(i) != ca->getSecondPosition(j) )
  {
    test2 = true;
    Nfix--;
  }
  if ( (!SiteTools::isConstant(s3in) || !SiteTools::isConstant(s3out)) && ca->getThirdPosition(i) != ca->getThirdPosition(j) )
  {
    test3 = true;
    Nfix--;
  }
  // Suppression of differences when not fixed
  vector<size_t> v(2);
  if (Nfix == 0)
  {
    v[0] = 0;
    v[1] = 0;
    return v;
  }
  if (Nfix < Ntot)
  {
    if (Na == 0)
      Ns = Nfix;
    if (Ns == 0)
      Na = Nfix;
    else
    {
      if (Ntot == 3)
      {
        if (Nfix == 1)
        {
          if (test1 && test2)
          {
            Na = 0; Ns = 1;
          }
          if (test1 && test3)
          {
            Na = 1; Ns = 0;
          }
          if (test2 && test3)
          {
            Na--; Ns--;
          }
        }
      }
      if (Nfix == 2)
      {
        if (test1)
        {
          Na = 1; Ns = 1;
        }
        if (test2)
          Na--;
        if (test3)
          Ns--;
      }
    }
    if (Ntot == 2)
    {
      if (test1)
      {
        if (ca->getSecondPosition(i) == ca->getSecondPosition(j))
          Na--;
        else
          Ns--;
      }
      if (test2)
        Na--;
      if (test3)
        Ns--;
    }
  }
  v[0] = Ns;
  v[1] = Na;
  return v;
}

/******************************************************************************/

bool CodonSiteTools::isFourFoldDegenerated(const Site& site, const GeneticCode& gc)
{
  if (!SiteTools::isConstant(site, true))
  {
    /** If non-synonymous mutation **/
    if (!(CodonSiteTools::isSynonymousPolymorphic(site, gc)))
      return false;

    for (size_t i = 0; i < site.size(); i++)
    {
      if (!(gc.isFourFoldDegenerated(site.getValue(i))))
      {
        return false;
      }
    }
  }
  else
  {
    for (size_t i = 0; i < site.size(); i++)
    {
      if (!(gc.isFourFoldDegenerated(site.getValue(i))))
      {
        return false;
      }
    }
  }
  return true;
}

/******************************************************************************/
// EGGLIB STUFF
/* ****************************************************************** */


    void BppDiversity::computeSingle(bpp::AlignedSequenceContainer& align) {

                      v_Ss = bpp::SequenceStatistics::polymorphicSiteNumber(align);
                    v_Sinf = bpp::SequenceStatistics::parsimonyInformativeSiteNumber(align);
                    v_Ssin = bpp::SequenceStatistics::countSingleton(align);
                     v_eta = bpp::SequenceStatistics::totNumberMutations(align);
                      v_He = bpp::SequenceStatistics::heterozygosity(align);
                     v_He2 = bpp::SequenceStatistics::squaredHeterozygosity(align);
                      //v_GC = bpp::SequenceStatistics::gcContent(align);
                      v_tW = bpp::SequenceStatistics::watterson75(align);
                     v_T83 = bpp::SequenceStatistics::tajima83(align);
                       v_K = bpp::SequenceStatistics::DVK(align);
                       v_H = bpp::SequenceStatistics::DVH(align);
                      v_Ti = bpp::SequenceStatistics::getNumberOfTransitions(align);
                      v_Tv = bpp::SequenceStatistics::getNumberOfTransversions(align);
   if (v_Tv!=0)     v_TiTv = bpp::SequenceStatistics::getTransitionsTransversionsRatio(align);
   else             v_TiTv = 0.;
   if (v_Ss!=0)        v_D = bpp::SequenceStatistics::tajimaDSS(align);
   else                v_D = 0.;
   if (v_eta!=0)    v_Deta = bpp::SequenceStatistics::tajimaDTNM(align);
   else             v_Deta = 0.;
   if (v_eta!=0) v_Dflstar = bpp::SequenceStatistics::fuliDstar(align);
   else          v_Dflstar = 0.;
   if (v_eta!=0)   v_Fstar = bpp::SequenceStatistics::fuliFstar(align);
   else            v_Fstar = 0.;
                    v_rhoH = bpp::SequenceStatistics::hudson87(align);

        if (b_coding) {
            v_ncodon1mut = bpp::SequenceStatistics::monoSitePolymorphicCodonNumber(align);
                 v_nstop = bpp::SequenceStatistics::stopCodonSiteNumber(align);
                  v_nsyn = bpp::SequenceStatistics::synonymousPolymorphicCodonNumber(align, *p_code);
                   v_tWS = bpp::SequenceStatistics::watterson75Synonymous(align, *p_code);
                  v_tWNS = bpp::SequenceStatistics::watterson75NonSynonymous(align, *p_code);
                   v_PiS = bpp::SequenceStatistics::piSynonymous(align, *p_code);
                  v_PiNS = bpp::SequenceStatistics::piNonSynonymous(align, *p_code);
                v_Ssites = bpp::SequenceStatistics::meanSynonymousSitesNumber(align, *p_code);
               v_NSsites = bpp::SequenceStatistics::meanNonSynonymousSitesNumber(align, *p_code);
                    v_SS = bpp::SequenceStatistics::synonymousSubstitutionsNumber(align, *p_code);
                   v_SNS = bpp::SequenceStatistics::nonSynonymousSubstitutionsNumber(align, *p_code);
        }
    }


    void BppDiversity::computeDouble(bpp::AlignedSequenceContainer& ingroup, bpp::AlignedSequenceContainer& outgroup) {
            v_Sext = bpp::SequenceStatistics::totMutationsExternalBranchs(ingroup, outgroup);
    if (v_eta!=0) v_Dfl = bpp::SequenceStatistics::fuliD(ingroup, outgroup);
                else v_Dfl = 0.;
    if (v_eta!=0) v_F = bpp::SequenceStatistics::fuliF(ingroup, outgroup);
                else v_F = 0.;

            if (b_coding) {
                v_MK = bpp::SequenceStatistics::MKtable(ingroup, outgroup, *p_code);
                v_NI = bpp::SequenceStatistics::neutralityIndex(ingroup, outgroup, *p_code);
            }
    }


/* ****************************************************************** */

    bool BppDiversity::hasOutgroup() const {
        return b_outgroup;
    }

    unsigned int BppDiversity::S() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Ss;           
    }
    
    unsigned int BppDiversity::Sinf() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Sinf;
    }

    unsigned int BppDiversity::Ssin() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Ssin;         
    }

    unsigned int BppDiversity::eta() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_eta;          
    }

    unsigned int BppDiversity::Sext() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        return v_Sext;         
    }

    double BppDiversity::He() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_He;           
    }

    double BppDiversity::He2() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_He2;          
    }

//    double BppDiversity::GC() const { 
//        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
//        return v_GC;           
//    }

    double BppDiversity::tW() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_tW;           
    }

    double BppDiversity::T83() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_T83;          
    }

    unsigned int BppDiversity::K() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_K;            
    }

    double BppDiversity::H() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_H;            
    }

    unsigned int BppDiversity::Ti() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Ti;           
    }

    unsigned int BppDiversity::Tv() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Tv;           
    }

    double BppDiversity::TiTv() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_TiTv;         
    }

    unsigned int BppDiversity::nstop() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_nstop;        
    }

    unsigned int BppDiversity::ncodon1mut() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_ncodon1mut;   
    }

    unsigned int BppDiversity::nsyn() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_nsyn;         
    }

    double BppDiversity::tWS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_tWS;          
    }

    double BppDiversity::tWNS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_tWNS;         
    }

    double BppDiversity::PiS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_PiS;          
    }

    double BppDiversity::PiNS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_PiNS;         
    }

    double BppDiversity::Ssites() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_Ssites;       
    }

    double BppDiversity::NSsites() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_NSsites;      
    }

    unsigned int BppDiversity::SS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_SS;           
    }

    unsigned int BppDiversity::SNS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_SNS;          
    }

    const std::vector<unsigned int>& BppDiversity::MK() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
       return v_MK;           
    }

    double BppDiversity::NI() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        return v_NI;           
    }

    double BppDiversity::D() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_D;            
    }

    double BppDiversity::Deta() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Deta;         
    }

    double BppDiversity::Dfl() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        return v_Dfl;          
    }

    double BppDiversity::Dflstar() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Dflstar;      
    }

    double BppDiversity::F() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        return v_F;            
    }

    double BppDiversity::Fstar() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Fstar;        
    }

    double BppDiversity::rhoH() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_rhoH;         
    }

}

/*

 DOCSTRING FOR egglib.Align.polymorphismBPP
 
 
 Computes diversity statistics using tools provided through the
 Bio++ libraries. Note that attempting to call this method from
 an EggLib module compile without Bio++ support will result in a
 RuntimeError.
 
 Arguments:
 
 * *dataType*: 1 for DNA, 2 for RNA, 3 for protein sequences,
 4 for standard codons, 5 for vertebrate mitochondrial
 codons, 6 for invertebrate mitochondrial codons and 7 for
 echinoderm mitochondrial codons.
 
 The method returns a dictionary containing the diversity
 statistics. Some keys will be computed only in the presence of
 an outgroup, or if sequences were specified as coding or
 depending on the value of other statistics (otherwise, they will
 be ``None``).
 
 The following statistics are always computed:
 
 * ``S``: Number of polymorphic sites.
 * ``Sinf``: Number of parsimony informative sites.
 * ``Ssin``: Number of singleton sites.
 * ``eta``: Minimal number of mutations.
 * ``thetaW``: Theta estimator (Watterson *Theor. Popul.
 Biol.* **7**:256-276, 1975).
 * ``T83``: Theta estimator (Tajima *Genetics* **105**:437-460, 1983)
 * ``He``: Heterozygosity.
 * ``Ti``: Number of transitions.
 * ``Tv``: Number of transversions.
 * ``K``: Number of haplotypes.
 * ``H``: Haplotypic diversity.
 * ``rhoH``: Hudson's estimator of rho (*Genet. Res.*
 **50**:245-250, 1987).
 
 The following statistic is computed only if ``Tv`` > 0:
 
 * ``TiTv``: Transition/transversion ratio.
 
 The following statistic is computed only if ``S`` > 0:
 
 * ``D``: Tajima statistic (*Genetics* **123**:585-595, 1989).
 
 The following statistics are computed only if ``eta`` > 0:
 
 * ``Deta``: Tajima's D computed with ``eta`` instead of ``S``.
 * ``Dflstar``: Fu and Li's D* (without outgroup; *Genetics* **133**:693-709).
 * ``Fstar``: Fu and Li's F* (without ougroup; *Genetics* **133**:693-709).
 
 The following statistic is computed only if an outgroup is found:
 
 * ``Sext``: Mutations on external branches.
 
 The following statistics are computed only if an outgroup is
 found and ``eta`` > 0:
 
 * ``Dfl``: Fu and Li's D (*Genetics* **133**:693-709).
 * ``F``: Fu and Li's F (*Genetics* **133**:693-709).
 
 The following statistics are computed only if sequences are
 coding *dataType* = 4-7:
 
 * ``ncodon1mut``: Number of codon sites with exactly one mutation.
 * ``NSsites``: Average number of non-synonymous sites.
 * ``nstop``: Number of codon sites with a stop codon.
 * ``nsyn``: Number of codon sites with a synonymous change.
 * ``PiNS``: Nucleotide diversity computed on non-synonymous sites.
 * ``PiS``: Nucleotide diversity computed on synonymous sites.
 * ``SNS``: Number of non-synonymous polymorphic sites.
 * ``SS``: Number of synonymous polymorphic sites.
 * ``Ssites``: Number of synonymous sites.
 * ``tWNS``: Watterson's theta computed on non-synonymous sites.
 * ``tWS``: Watterson's theta computed on synonymous sites.
 
 The following statistics are computed only if sequences are
 coding *dataType* = 4-7 and an outgroup is found:
 
 * ``MK``: McDonald-Kreitman test table (*Nature*
 **351**:652-654, 1991).
 * ``NI``: Neutrality index (Rand and Kann *Mol. Biol. Evol.*
 **13**:735-748).
 
 The returned dictionary also contains a nest dictionary ``options``
 which feedbacks the values used at function call.
 
 .. versionchanged:: 2.0.2
 The following statistics are now computed only if ``S`` > 0:
 ``D``, ``Deta``, ``Dflstar``, ``Fstar``, ``Dfl``, ``F``.
 
 .. versionchanged:: 2.1.0
 The statistics not computed are now exported and set to
 ``None``.
 
*/