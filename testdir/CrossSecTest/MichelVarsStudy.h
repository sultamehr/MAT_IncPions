for(size_t whichMichel = 0; whichMichel < myevent.m_nmichels.size(); ++whichMichel)
            {
              (*var->m_bestPionByGENIELabel)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValue(*universe, whichMichel), universe->GetWeight());
            }
