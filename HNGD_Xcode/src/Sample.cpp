#include "Sample.hpp"
#include "PhysicsConstants.h"
#include <iostream>

Sample :: Sample(int nbCells, double bias, double tssp0, double Qp, double tssd0, double Qd) :

    _nbCells(nbCells),
    _bias   (bias),

    _tssp0  (tssp0),
    _Qp     (Qp),

    _tssd0  (tssd0),
    _Qd     (Qd)
{
    _position       = vector<double>(_nbCells) ;
    _temperature    = vector<double>(_nbCells) ;
    _solutionContent= vector<double>(_nbCells) ;
    _hydrideContent = vector<double>(_nbCells) ;
    _totalContent   = vector<double>(_nbCells) ;
    _tssd           = vector<double>(_nbCells) ;
    _tssp           = vector<double>(_nbCells) ;
}

// Compute the equilibrium for the initial conditions
void Sample :: computeEquilibrium()
{
    for(int k=0; k<_nbCells; k++)
    {
        _solutionContent[k] = min(_totalContent[k], _tssd[k]) ;
        _hydrideContent[k] = _totalContent[k] - _solutionContent[k] ;
    }
}


// Solubilities computation
void Sample :: computeTSS()
{
    for(int k=0; k<_nbCells; k++)
    {
        _tssd[k] = _tssd0 * exp(-_Qd / (R * _temperature[k])) ;
        _tssp[k] = _tssp0 * exp(-_Qp / (R * _temperature[k])) ;
    }
    
}

// Domain definition
void Sample :: computeLocations(double x0, double xEnd, int _geometry)
{
	if (_geometry > 0){ // Polar, want to not create an Xend point since 2pi = 0 radians
		double xEnd = 2*M_PI;
	    const double initialLenght = 2*M_PI/_nbCells;

	    _position[0] = x0 ;
	    for (int k=1; k<_nbCells; k++)
	    	_position[k] = _position[k-1] + initialLenght;

	    _position[_nbCells] = xEnd - initialLenght;
    }

	else { // Linear
	    double  sum = 1. + _bias ;
	    for(int k=0; k<_nbCells-3; k++)
	        sum = 1. + _bias*sum ;

	    const double initialLenght = (xEnd - x0)/sum ;

	    _position[0] = x0 ;
	    _position[1] = x0 + initialLenght ;

	    for(int k=2; k<_nbCells-1; k++)
	        _position[k] = _position[k-1] + _bias*(_position[k-1] - _position[k-2]) ;

	    _position[_nbCells-1] = xEnd  ;}
}

// Interpolation
void Sample :: spatialeInterpolation(vector<double>& refX, vector<double>& refY, vector<double>& vectorY)
{
    int k = 1 ;
    for(int i=0; i<_position.size(); i++) {
        if(_position[i] > refX[k])
            k++ ;
        vectorY[i] = refY[k-1] + (refY[k] - refY[k-1]) * (_position[i] - refX[k-1]) / (refX[k] - refX[k-1]);
    }
}

void Sample :: polarInterpolation(vector<double>& refX, vector<double>& refY, vector<double>& vectorY)
{
	vector <double> localposition;
	localposition = _position;
	int k = 1;
	double initialLenght = 2*M_PI / localposition.size() ;
	localposition.push_back(localposition[localposition.size()]+initialLenght);

	for(int i=0; i<localposition.size(); i++) {
		if(localposition[i] > refX[k])
			k++ ;
		vectorY[i] = refY[k-1] + (refY[k] - refY[k-1]) * (localposition[i] - refX[k-1]) / (refX[k] - refX[k-1]);
	}
}
