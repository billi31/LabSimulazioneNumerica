//
//  random.hpp
//  lab2
//
//  Created by Administrator on 26/03/21.
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
	double cos_IS();
	void Walk_d(double &, double &, double &);
	void Walk_c(double &, double &, double &);
	
};

#endif // __Random__
