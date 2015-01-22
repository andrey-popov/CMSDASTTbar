#include <Reader.hpp>

#include <stdexcept>
#include <sstream>
#include <algorithm>


using namespace std;


// A static data member
unsigned const Reader::maxSize;


Reader::Reader(shared_ptr<TFile> &srcFile_, list<string> const &treeNames_, bool isMC_ /*= true*/):
    srcFile(srcFile_), treeNames(treeNames_), curTreeNameIt(treeNames.begin()), isMC(isMC_),
    curSystType(SystType::Nominal), curSystDirection(SystDirection::Up),
    applyBTagReweighting(true)
{
    // Make sure the source file is a valid one
    if (not srcFile or srcFile->IsZombie())
        throw runtime_error("The source file does not exist or is corrupted.");
    
    
    // Get the first tree
    GetTree(*curTreeNameIt);
}


Reader::Reader(shared_ptr<TFile> &srcFile_, string const &treeName, bool isMC_ /*= true*/):
    Reader(srcFile_, list<string>{treeName}, isMC_)
{}


bool Reader::ReadNextEvent()
{
    // Check if there are events left in the current source tree
    if (curEntry == nEntries)  // no more events in the current tree
    {
        ++curTreeNameIt;
        
        if (curTreeNameIt == treeNames.end())  // no more source trees
            return false;
        
        GetTree(*curTreeNameIt);
    }
    
    
    // Either there were events in the current source file or a new file has been opened
    curTree->GetEntry(curEntry);
    ++curEntry;
    
    
    // Copy properies of objects in the event from read buffers
    leptons.clear();
    
    for (int i = 0; i < lepSize; ++i)
        leptons.emplace_back(lepFlavour[i], lepPt[i], lepEta[i], lepPhi[i], lepIso[i]);
    
    jets.clear();
    
    for (int i = 0; i < jetSize; ++i)
        jets.emplace_back(jetPt[i], jetEta[i], jetPhi[i], jetBTag[i], jetFlavour[i]);
    
    met.Set(metPt, metPhi);
    
    if (isMC)
    {
        jetsJECUp.clear();
        jetsJECDown.clear();
        
        for (int i = 0; i < jetJECUpSize; ++i)
            jetsJECUp.emplace_back(jetJECUpPt[i], jetJECUpEta[i], jetJECUpPhi[i],
             jetJECUpBTag[i], jetJECUpFlavour[i]);
        
        for (int i = 0; i < jetJECDownSize; ++i)
            jetsJECDown.emplace_back(jetJECDownPt[i], jetJECDownEta[i], jetJECDownPhi[i],
             jetJECDownBTag[i], jetJECDownFlavour[i]);
        
        metJECUp.Set(metJECUpPt, metJECUpPhi);
        metJECDown.Set(metJECDownPt, metJECDownPhi);
    }
    
    
    // Make sure vector of leptons and jets are ordered in pt
    sort(leptons.rbegin(), leptons.rend());
    sort(jets.rbegin(), jets.rend());
    
    if (isMC)
    {
        sort(jetsJECUp.rbegin(), jetsJECUp.rend());
        sort(jetsJECDown.rbegin(), jetsJECDown.rend());
    }
    
    
    // Indicate that the stored event weight is no longer up-to-date
    weightCached = false;
        
    
    return true;
}


void Reader::Rewind() noexcept
{
    curTreeNameIt = treeNames.begin();
    
    // Delete the current tree. It must be done during the rewind because if there is only single
    //tree in the sample, the GetTree will try to reset curTree to the same pointer, and it will
    //lead to a segfault
    curTree.reset();
    
    GetTree(*curTreeNameIt);
}


void Reader::SetSystematics(SystType systType, SystDirection systDirection)
{
    // Update information about requested systematics
    curSystType = systType;
    curSystDirection = systDirection;
    
    
    // If the type is Nominal (i.e. no variation), only direction Up is allowed
    if (curSystType == SystType::Nominal)
        curSystDirection = SystDirection::Up;
    
    
    // The stored weight might not be up-to-date anymore since it might be affected by the
    //systematics
    weightCached = false;
}


vector<Lepton> const &Reader::GetLeptons() const noexcept
{
    return leptons;
}


vector<Jet> const &Reader::GetJets() const noexcept
{
    if (isMC and curSystType == SystType::JEC)
    {
        if (curSystDirection == SystDirection::Up)
            return jetsJECUp;
        else
            return jetsJECDown;
    }
    else
        return jets;
}


MET const &Reader::GetMET() const noexcept
{
    if (isMC and curSystType == SystType::JEC)
    {
        if (curSystDirection == SystDirection::Up)
            return metJECUp;
        else
            return metJECDown;
    }
    else
        return met;
}


double Reader::GetWeight() noexcept
{
    // If the current sample is data, the answer is trivial
    if (not isMC)
        return 1.;
    
    
    // Check if the weight is up-to-date
    if (weightCached)
        return weight;
    
    
    // Recalculate the weight. Note that if the workflow reaches this point, the current sample is
    //simulation
    
    // Raw weights stored in the trees inlcude effects of pile-up, lepton scale factors, and
    //normalisation for the cross section and integrated luminosity
    weight = rawWeight;
    
    
    // Reweighting for the b-tagging scale factors
    if (applyBTagReweighting)
        for (auto const &j: jets)
        {
            double const perJetBTagWeight =
             csvReweighter.CalculateJetWeight(j, curSystType, curSystDirection);
            
            if (perJetBTagWeight != 0.)
                weight *= perJetBTagWeight;
        }
    
    
    // The weight is now cached
    weightCached = true;
    
    
    // Return the weight
    return weight;
}


unsigned Reader::GetNumPV() const noexcept
{
    return nPV;
}


void Reader::SwitchBTagReweighting(bool on /*= true*/)
{
    applyBTagReweighting = on;
}


void Reader::GetTree(string const &name)
{
    // Get the tree from the source file
    //curTree.reset(dynamic_cast<TTree *>(srcFile->Get(name.c_str())));
    curTree.reset(dynamic_cast<TTree *>(srcFile->Get(name.c_str())));
    
    
    // Make sure the tree exists
    if (not curTree)
    {
        ostringstream ost;
        ost << "Cannot find tree \"" << name << "\" in file \"" << srcFile->GetTitle() << "\".";
        throw runtime_error(ost.str());
    }
    
    
    // Set event counters
    nEntries = curTree->GetEntries();
    curEntry = 0;
    
    
    // Set buffers to read the tree
    curTree->SetBranchAddress("nlepton", &lepSize);
    curTree->SetBranchAddress("lept_pt", lepPt);
    curTree->SetBranchAddress("lept_eta", lepEta);
    curTree->SetBranchAddress("lept_phi", lepPhi);
    curTree->SetBranchAddress("lept_iso", lepIso);
    curTree->SetBranchAddress("lept_flav", lepFlavour);
    
    curTree->SetBranchAddress("njets", &jetSize);
    curTree->SetBranchAddress("jet_pt", jetPt);
    curTree->SetBranchAddress("jet_eta", jetEta);
    curTree->SetBranchAddress("jet_phi", jetPhi);
    curTree->SetBranchAddress("jet_btagdiscri", jetBTag);
    curTree->SetBranchAddress("jet_flav", jetFlavour);
    
    curTree->SetBranchAddress("met_pt", &metPt);
    curTree->SetBranchAddress("met_phi", &metPhi);
    
    curTree->SetBranchAddress("nvertex", &nPV);
    
    if (isMC)
    {
        curTree->SetBranchAddress("jesup_njets", &jetJECUpSize);
        curTree->SetBranchAddress("jet_jesup_pt", jetJECUpPt);
        curTree->SetBranchAddress("jet_jesup_eta", jetJECUpEta);
        curTree->SetBranchAddress("jet_jesup_phi", jetJECUpPhi);
        curTree->SetBranchAddress("jet_jesup_btagdiscri", jetJECUpBTag);
        curTree->SetBranchAddress("jet_jesup_flav", jetJECUpFlavour);
        
        curTree->SetBranchAddress("jesdown_njets", &jetJECDownSize);
        curTree->SetBranchAddress("jet_jesdown_pt", jetJECDownPt);
        curTree->SetBranchAddress("jet_jesdown_eta", jetJECDownEta);
        curTree->SetBranchAddress("jet_jesdown_phi", jetJECDownPhi);
        curTree->SetBranchAddress("jet_jesdown_btagdiscri", jetJECDownBTag);
        curTree->SetBranchAddress("jet_jesdown_flav", jetJECDownFlavour);
        
        curTree->SetBranchAddress("met_jesup_pt", &metJECUpPt);
        curTree->SetBranchAddress("met_jesup_phi", &metJECUpPhi);
        
        curTree->SetBranchAddress("met_jesdown_pt", &metJECDownPt);
        curTree->SetBranchAddress("met_jesdown_phi", &metJECDownPhi);
        
        curTree->SetBranchAddress("evtweight", &rawWeight);
    }
    
    
    // Set the event weight for data (it will not be modified)
    weight = 1.;
}
