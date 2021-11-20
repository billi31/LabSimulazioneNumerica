//
//  random.hpp
//  lab1
//
//  Created by Administrator on 25/03/21.
//
//


#ifndef __Random__
#define __Random__

class Random {
	
private:
	int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;
	
protected:
	
public:
	// constructors
	Random();
	// destructor
	~Random();
	// methods
	void Initialize(Random &);
	void SetRandom(int * , int, int);
	void SaveSeed();
	double Rannyu(void);
	double Rannyu(double min, double max);
	double Gauss(double mean, double sigma);
	double exponential(double lambda);
	double CauchyLorentz(double mean, double Gamma);
	
};

#endif // __Random__
