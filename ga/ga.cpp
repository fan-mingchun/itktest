
/*
* 
* website: http://www.coolsoft-sd.com/
* contact: support@coolsoft-sd.com
*
*/

#include <iostream>

#include "Initialization.h"
#include "ChromosomeOperations.h"
#include "MultiValueChromosome.h"
#include "Population.h"
#include "StopCriterias.h"
#include "IncrementalAlgorithm.h"

using namespace std;

using namespace Algorithm;
using namespace Algorithm::StopCriterias;
using namespace Algorithm::SimpleAlgorithms;

bool initParaFlag = false;
class fFitness : public GaFitnessOperation
{
public:

	virtual float GACALL operator ()(const GaChromosome* chromosome) const
	{
		auto chrome = const_cast<GaChromosome*>(chromosome);
		vector<double>& vals = dynamic_cast< GaMVArithmeticChromosome<double>*>( chrome )->GetCode(0);
		/*double* initp = new double[10];
		initp[0] = 50;*/


		if(!initParaFlag)
		{
			//��ʼ������
			initParaFlag = true;
			int prinNum = 10;
			for(int index = 0;index < prinNum;++index)
			{
				vals[0] = 0;
			}
			//pos
			vals[prinNum++]=(13/91.0)*2-1;
			vals[prinNum++]=(57/188.0)*2-1;
			vals[prinNum++]=(20/19.0)*2-1;
			//angle
			vals[prinNum++]=0;
			vals[prinNum++]=0;
			vals[prinNum++]=0;

			vals[prinNum++]=(1-0.8)*5-1;
			

		}
		int pos[3];
		pos[0] = ((vals[10]+1)/2)*91;
		pos[1] = ((vals[11]+1)/2)*188;
		pos[2] = ((vals[12]+1)/2)*19;
		double angle[3];
		angle[0] = vals[13];
		angle[1] = vals[14];
		angle[2] = vals[15];

		double scale = (vals[16]+1)/5+0.8;

		return (float)((100-vals[0])*vals[0] );
	}

	virtual GaParameters* GACALL MakeParameters() const { return NULL; }

	virtual bool GACALL CheckParameters(const GaParameters& parameters) const { return true; }
};

class fObserver : public GaObserverAdapter
{
	virtual void GACALL NewBestChromosome(const GaChromosome& newChromosome, const GaAlgorithm& algorithm)
	{
		const vector<double>& vals = dynamic_cast<const GaMVArithmeticChromosome<double>&>( newChromosome ).GetCode();
		cout << "New chromosome found:\n";
		cout << "Fitness: " << newChromosome.GetFitness() << endl;
		cout << vals[0] << endl;
	}

	virtual void GACALL EvolutionStateChanged(GaAlgorithmState newState, const GaAlgorithm& algorithm)
	{
		if( newState == GAS_RUNNING )
			cout << "start\n";
		else if( newState == GAS_CRITERIA_STOPPED )
			cout << "end";
	}
};

int main()
{
	const int paraNum = 17;
	GaInitialize();

	GaChromosomeParams chromosomeParams( 0.03f, 1, true, 0.8f, 1 );

	GaValueIntervalBounds<double> valueInt( -1, 1 );
	GaValueIntervalBounds<double> invValueInt( -1, 1 );
	GaIntervalValueSet<double> valueSet( valueInt, invValueInt, GaGlobalRandomDoubleGenerator, false);

	fFitness fitnessOperation;
	GaChromosomeDomainBlock<double> configBlock( &valueSet,
		GaCrossoverCatalogue::Instance().GetEntryData( "GaMultiValueCrossover" ),
		GaMutationCatalogue::Instance().GetEntryData( "GaFlipMutation" ),
		&fitnessOperation, GaFitnessComparatorCatalogue::Instance().GetEntryData( "GaMaxFitnessComparator" ),
		&chromosomeParams );

	GaMVArithmeticChromosome<double> prototype( paraNum, &configBlock );

	GaPopulationConfiguration populationConfig;
	GaPopulationParameters populationParams( 30, false, true, false, 0, 0 );

	populationConfig.SetParameters( populationParams );
	populationConfig.SetSortComparator( &configBlock.GetFitnessComparator() );

	GaPopulation population( &prototype, &populationConfig );


	GaMultithreadingAlgorithmParams algorithmParams( 1 );
	Algorithm::SimpleAlgorithms::GaIncrementalAlgorithm algorithm( &population, algorithmParams );

	GaGenerationCriteriaParams criteriaParams( 100000 );
	algorithm.SetStopCriteria( GaStopCriteriaCatalogue::Instance().GetEntryData( "GaGenerationCriteria" ), &criteriaParams);

	fObserver observer;
	algorithm.SubscribeObserver( &observer );

	algorithm.StartSolving( false );

	algorithm.WaitForThreads();

	GaFinalize();

	return 0;
}
