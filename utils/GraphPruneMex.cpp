#include "mex.h"
#include "matrix.h"
#include "math.h"
#include <stdlib.h>
#include <cstring>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;


template <class Type>
class Matrix
{
public:
	typedef Type* iterator;

	Type* data;
	size_t rows;
	size_t cols;
	
	Matrix()
	{
		rows = 0;
		cols = 0;
		data = NULL;
	}
	Matrix( const Matrix<Type> & mat1 )
	{
		rows = mat1.rows;
		cols = mat1.cols;
		if ( size() > 0)
		{
			data = new Type [rows*cols];
			memcpy(data, mat1.data, sizeof(Type) * rows * cols);
		}
	}
	Matrix( const mxArray * array1 )
	{
		rows = mxGetM(array1);
		cols = mxGetN(array1);
		if (size() > 0)
		{
			double* ptr = mxGetPr(array1);
			data = new Type [rows*cols];
			for (int i = 0; i < size(); i++)
				data[i] = ptr[i];
		}
	}
	~Matrix()
	{
		if (size() > 0)
			delete [] data;
	}
	void Create(size_t s1, size_t s2)
	{
		rows = s1;
		cols = s2;
		if (size() > 0)
		{
			data = new Type [rows*cols];
			memset(data, 0, sizeof(Type) * rows * cols);
		}
	}
	Type * begin()
	{	
		return data;	
	}
	Type * end()
	{
		return data + rows * cols;	
	}
	size_t size()
	{
		return rows * cols;
	}
	Type & operator[] ( size_t index )
	{
		return data[index];
	}
	void Set(Type value)
	{
		for (int i = 0; i < rows * cols; i++)
			data[i] = value;
	}
	void Add(Type value)
	{
		for (int i = 0; i < rows * cols; i++)
			data[i] += value;
	}
	void Copy( Matrix<Type> & mat1, int start, int end )
	{
		for (int i = start; i < end; i++)
			data[i] = mat1[i];
	}

	template <class SType>
	void Copy( set<SType> & set1, size_t start = 0 )
	{
		typename set<SType>::iterator it = set1.begin();
		for (int i = start; i < start + set1.size(); i++, it++)
			data[i] = *it;
	}

	int Shrink(int maxID)
	{
		// Renumerate region ID
		int * mark = new int [maxID+1];
		memset(mark, 0, sizeof(int) * (maxID+1));
		for (int i = 0; i < rows; i++)
			mark[data[i]] = 1;

		int num_zero = 0;
		int * label_trans = new int [maxID + 1];
		for (int i = 0; i < maxID + 1; i++)
		{
			if (mark[i] == 0)
				num_zero++;
			else
				label_trans[i] = i - num_zero;
		}

		for (int i = 0; i < rows; i++)
			data[i] = label_trans[data[i]];
		
		delete [] label_trans;
		delete [] mark;
		
		return maxID + 1 - num_zero;
	}
};
#define EPS 1e-10

bool pairless( const pair<int, int>& i, const pair<int, int>& j )
{
	return ( i.first < j.first );
}

// Blist and supp should be in ascending order
void PruneBlocktmp( Matrix<double> & cv,  Matrix<int> & Bm,  Matrix<int> & BC,  double lamada,
					double cl,        Matrix<int> & supp,     int RegionNum,        Matrix<int> & Blist,      Matrix<int> & RegionIndex, 
					double & cl_out,  Matrix<int> & supp_out, int & RegionNum_out,  Matrix<int> & Blist_out,  Matrix<int> & RegionIndex_out,
					Matrix<int> & new_supp)
{
	// Remove the blocks in Blist 
    Matrix<int> curr_possBlist;
	curr_possBlist.Create(Bm.rows, 1);
	curr_possBlist.Set(1);
	for (int i = 0; i < Blist.size(); i++)
		curr_possBlist[Blist[i]] = -1;

	// Remove the blocks if all elements in this block are in supp already
	Matrix<int> newindexnum;
	newindexnum.Create(Bm.rows, 1);
	for (int i = 0; i < Bm.rows; i++)
	{
		if (curr_possBlist[i] < 0)
			continue;
		int num_no_supp = 0;
		for (int j = 0; j < Bm.cols; j++)
		{
			if ( Bm[j*Bm.rows+i] >= 0 && !binary_search(supp.begin(), supp.end(), Bm[j*Bm.rows+i]) )	// not in supp
				num_no_supp ++;
		}
		if (num_no_supp > 0)
			newindexnum[i] = num_no_supp;
		else
			curr_possBlist[i] = -1;
	}
		
	// Compute the increase of region number: delta_cl
	Matrix<double> delta_cl;
	delta_cl.Create(BC.rows, 1);
	delta_cl.Set(0);
	for (int i = 0; i < BC.rows; i++)
	{	
		if (curr_possBlist[i] < 0)
			continue;
		set<int> regions;
		for (int j = 0; j < BC.cols; j++)
		{
			pair<int*, int*> pos = equal_range( Blist.begin(), Blist.end(), BC[j*BC.rows+i]);
			if (  pos.first != pos.second )
			{
				int region = RegionIndex[pos.first - Blist.begin()];
				regions.insert(region);
			}
		}
		delta_cl[i] = 1 - (int) regions.size();
	}
	
	// delta_cl=delta_cl1*log2(m)+lamada*newindexnum; 
	double k = log(double(Bm.rows))/log(2.0);
	for (int i = 0; i < BC.rows; i++)
	{
		if (curr_possBlist[i] < 0)
			continue;
		delta_cl[i] = lamada * delta_cl[i] * k + newindexnum[i];
	}

	// Compute error
	Matrix<double> lqerror;
	lqerror.Create(Bm.rows, 1);
	lqerror.Set(0);
	for (int i = 0; i < Bm.rows; i++)
	{
		if (curr_possBlist[i] < 0)
			continue;
		for (int j = 0; j < Bm.cols; j++)
		{
			if ( Bm[j*Bm.rows+i] < 0 || binary_search( supp.begin(), supp.end(), Bm[j*Bm.rows+i]) )		// in supp
				continue;
			double err = cv[Bm[j*Bm.rows+i]];
			lqerror[i] += err * err;
		}
	}
	
	// Select new block
	int Select_bind;
	double min_delta_cl = *min_element(delta_cl.begin(), delta_cl.end());
	if (min_delta_cl >= 0)
	{
		for (int i = 0; i < Bm.rows; i++)
		{
			if (curr_possBlist[i] < 0)
				lqerror[i] = 0;
			else
				lqerror[i] /= delta_cl[i];
		}
		Select_bind = max_element(lqerror.begin(), lqerror.end()) - lqerror.begin();
	}
	else
	{
		double max_lqerror = 0;
		int max_lqerror_id = -1;
		for (int i = 0; i < BC.rows; i++)
			if (delta_cl[i] <= min_delta_cl + EPS)
			{
				if (max_lqerror < lqerror[i])
				{
					max_lqerror = lqerror[i];
					max_lqerror_id = i;
				}
			}
		Select_bind = max_lqerror_id;
	}
		
	///////////////////////////////////////////////////////////////
	cl_out = cl + delta_cl[Select_bind];
	
	// Create list of newly added support nodes
	set<int> new_supp_set;
	for (int j = 0; j < Bm.cols; j++)
	{
		int bid = Bm[ j * Bm.rows + Select_bind];
		if ( bid >= 0 && !binary_search( supp.begin(), supp.end(), bid) )
			new_supp_set.insert(bid);
	}
	new_supp.Create( new_supp_set.size(), 1 );
	new_supp.Copy( new_supp_set );
	supp_out.Create( supp.size() + new_supp_set.size() , 1 );
	supp_out.Copy(supp, 0, supp.size());
	supp_out.Copy(new_supp_set, supp.size());
	sort(supp_out.begin(), supp_out.end());

	// Create new Blist
	Blist_out.Create(Blist.rows+1, 1);
	Blist_out.Copy(Blist, 0, Blist.rows);
	Blist_out[Blist.rows] = Select_bind;
	
	///////////////////////////////////////////////////////////////
	// Merge regions
	int * label_trans = new int [RegionNum + 1];
	for (int i = 0; i < RegionNum + 1; i++)
		label_trans[i] = i;
	for (int j = 0; j < BC.cols; j++)
	{
		int bid = BC[ j * BC.rows + Select_bind];
		pair<int*, int*> pos = equal_range( Blist.begin(), Blist.end(), bid);
		if (  pos.first != pos.second )
		{
			int region = RegionIndex[pos.first - Blist.begin()];
			label_trans[region] = RegionNum;
		}
	}

	RegionIndex_out.Create(RegionIndex.rows + 1, 1);
	for ( int i = 0; i < RegionIndex.rows; i++ )
		RegionIndex_out[i] = label_trans[RegionIndex[i]];
	RegionIndex_out[RegionIndex.rows] = label_trans[RegionNum];

	RegionNum_out = RegionIndex_out.Shrink(RegionNum);
	delete [] label_trans;

	///////////////////////////////////////////////////////////////
	// Sort Blist_out and corresponding RegionIndex_out
	vector<pair<int, int> > Blist_RegionIndex;
	Blist_RegionIndex.resize(Blist_out.size());
	for (int i = 0; i< Blist_out.size(); i++)
		Blist_RegionIndex[i] = make_pair( Blist_out[i], RegionIndex_out[i] );
	
	sort(Blist_RegionIndex.begin(), Blist_RegionIndex.end(), pairless);
	for (int i = 0; i< Blist_out.size(); i++)
	{
		Blist_out[i] = Blist_RegionIndex[i].first;
		RegionIndex_out[i] = Blist_RegionIndex[i].second;
	}

}

template <class Type>
mxArray* CreateFromMat(Matrix<Type> & mat)
{
	mxArray* array1 = mxCreateDoubleMatrix(mat.rows, mat.cols, mxREAL);
	double* ptr = mxGetPr(array1);
	for (int i = 0; i < mat.size(); i++)
		ptr[i] = mat[i];
	return array1;
}

mxArray* CreateFromScalar(double scalar)
{
	mxArray* array1 = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* ptr = mxGetPr(array1);
	ptr[0] = scalar;
	return array1;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs!=9) 
		mexErrMsgTxt("9 inputs required: cv, Bm, BC, lamada, supp, Blist, cl0, RegionNum, RegionIndex");
	if(nlhs!=6) 
		mexErrMsgTxt("Two output required.");

	Matrix<double> cv(prhs[0]);
	Matrix<int>    Bm(prhs[1]);
	Matrix<int>    BC(prhs[2]);
	if ( BC.rows != Bm.rows )
		mexErrMsgTxt("Row num of Bm and Bc shoulbe be same");
	double lamada = mxGetScalar(prhs[3]);

	Matrix<int>    supp(prhs[4]);
	Matrix<int>    Blist(prhs[5]);
	double cl        = mxGetScalar(prhs[6]);
	int    RegionNum = mxGetScalar(prhs[7]);
	Matrix<int>    RegionIndex(prhs[8]);

	Bm.Add(-1);
	BC.Add(-1);
	supp.Add(-1);
	Blist.Add(-1);
	RegionIndex.Add(-1);

	double cl_out;
	int RegionNum_out;
	Matrix<int>    supp_out;
	Matrix<int>    Blist_out;
	Matrix<int>    RegionIndex_out;
	Matrix<int>    newindex;
	
	PruneBlocktmp( cv, Bm, BC, lamada,
					cl,      supp,     RegionNum,     Blist,     RegionIndex, 
					cl_out,  supp_out, RegionNum_out, Blist_out, RegionIndex_out,
					newindex );

	supp_out.Add(1);
	Blist_out.Add(1);
	RegionIndex_out.Add(1);
	newindex.Add(1);

	plhs[0] = CreateFromMat(supp_out);
	plhs[1] = CreateFromMat(Blist_out);
	plhs[2] = CreateFromScalar(cl_out);
	plhs[3] = CreateFromScalar(RegionNum_out);
	plhs[4] = CreateFromMat(RegionIndex_out);
	plhs[5] = CreateFromMat(newindex);
}
