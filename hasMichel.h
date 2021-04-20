//#include "CVUniverse.h"
//#include "Michel.h"

//PlotUtils includes
#include "PlotUtils/Cut.h"

template <class UNIVERSE, class EVENT>
class hasMichel: public PlotUtils::Cut<UNIVERSE, EVENT>
{
 public:
    hasMichel(): PlotUtils::Cut<UNIVERSE, EVENT>("Event Has Michel ")
    {
    }

  private:
    using Michel = typename std::remove_reference<decltype(*std::declval<EVENT>().m_nmichels.front())>::type;

    bool checkCut(const UNIVERSE& univ, EVENT& evt) const
    {
      int nmichels = univ.GetNMichels();
      evt.m_bestdist = 9999.;
      for (int i = 0; i < nmichels; ++i)
      {
        Michel* current_michel = new Michel(univ, i);
        if (current_michel->is_fitted != 1) continue;
         current_michel->DoesMichelMatchVtx(univ, current_michel);
         current_michel->DoesMichelMatchClus(univ, current_michel);

         current_michel->GetBestMatch(current_michel);


         double dist = current_michel->Best3Ddist;
           if (dist <= evt.m_bestdist) {
           evt.m_bestdist = dist;
           evt.m_idx = i;
           //std::cout << "Printing 2D distanecs" << current_michel->best_dist2D[0] << ", " <<  current_michel->best_dist2D[1] << ", "  <<  current_michel->best_dist2D[2] << std::endl;
           evt.m_best_XZ = current_michel->best_XZ;
           evt.m_best_UZ = current_michel->best_UZ;
           evt.m_best_VZ = current_michel->best_VZ;

	   }

         evt.m_nmichels.push_back(current_michel);
        
       }


        return !evt.m_nmichels.empty();
    }

};
