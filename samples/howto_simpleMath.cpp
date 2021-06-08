#include "simpleMath.h"

#include <iostream>
using namespace std;


int howTo_Matrix()
{
	cout << "\n*** Fill matrix with Parenthesis:" << endl;
	cout << "\t >> Matrix m(2, 2);" << endl;
	cout << "\t >> m(0,0) = 1;" << endl;

	Matrix m(2, 2);
	m(0, 0) = 1;
	m(0, 1) = 2;
	m(1, 0) = 3;
	m(1, 1) = 4;
	cout << m << endl;

	Matrix m2(2, 2);
	m2(0, 0) = 10;
	m2(0, 1) = 10;
	m2(1, 0) = 10;
	m2(1, 1) = 10;

	cout << "\n*** Print matrix: \n" << endl;
	cout << "\t >> cout << m << endl; \n" << endl;

	system("pause");

	cout << "\n*** Copy matrix (not the pointer): \n" << endl;
	cout << "\t >> Matrix m2copy(2,2); \n" << endl;
	cout << "\t >> Matrix m_copied = m2copy; \n" << endl;

	Matrix m2copy(2, 2);
	m2copy(0, 0) = 1;
	m2copy(0, 1) = 2;
	m2copy(1, 0) = 3;
	m2copy(1, 1) = 4;

	Matrix m_copied = m2copy;

	cout << m_copied << endl;

	cout << "\n*** Copy matrix (not the pointer): \n" << endl;
	cout << "\t >> Matrix *m2copy_p = new Matrix(2,2); \n" << endl;
	cout << "\t >> Matrix *m_copied_p = new Matrix(m2copy_p); \n" << endl;

	Matrix *m2copy_p = new Matrix(2, 2);
	(*m2copy_p)(0, 0) = 1;
	(*m2copy_p)(0, 1) = 2;
	(*m2copy_p)(1, 0) = 3;
	(*m2copy_p)(1, 1) = 4;

	Matrix *m_copied_p = new Matrix(m2copy_p);

	cout << m_copied_p << endl;

	delete m2copy_p;
	
	cout << m_copied_p << endl;	


	system("pause");

	cout << "\n*** Fill matrix with the constructor: \n" << endl;
	cout << "\t >> double datos3[4] = { 1, 2, 3, 4 };" << endl;
	cout << "\t >> Matrix m3(2, 2, datos); \n" << endl;

	double datos3[4] = { 1, 2, 3, 4 };
	Matrix m3(2, 2, datos3);
	cout << m3 << endl;
	m3(1, 1) = 10;
	cout << m3 << endl;

	system("pause");

	cout << "\n*** Fill square matrix with the constructor: \n" << endl;
	cout << "\t >> double datos[4] = { 1, 2, 3 };" << endl;
	cout << "\t >> Matrix m(2, datos); \n" << endl;

	Matrix m4(2, datos3);
	cout << m4 << endl;

	system("pause");

	cout << "\n*** Standard matrices: \n" << endl;
	cout << "\t >> m.Zeros();" << endl;
	Matrix m5(3, 3);
	m5.Zeros();
	cout << m5 << endl;
	cout << "\t >> m.Identity(); \n" << endl;
	m5.Identity();
	cout << m5 << endl;

	cout << "\t >> m.Ones(); \n" << endl;
	m5.Ones();
	cout << m5 << endl;

	cout << "\t >> m.Rand(); // TODO \n" << endl;
	//m5.Rand();
	//cout << m5 << endl;

	system("pause");

	cout << "\n*** Create Matrix equal to existing one: \n" << endl;
	cout << "\t >> Matrix m2 = m1;" << endl;
	Matrix m6 = m5;
	cout << m6 << endl;

	cout << "\t >> Matrix m2(m1);" << endl;
	Matrix m10(m5);
	cout << m10 << endl;

	system("pause");

	Matrix mout;
	cout << "\n*** Matrix Transpose: \n" << endl;
	cout << "\t >> m.Transpose(&mout); // The solution is stored in mout" << endl;
	cout << "Transpose(" << m << ")" << endl;
	m.Transpose(&mout);
	cout << mout << endl;
	cout << "\t >> m.Transpose(); // The solution is stored in m10" << endl;
	m.Transpose();
	cout << m << endl;


	system("pause");

	cout << "\n*** Matrix MULTIPLICATION: \n" << endl;
	cout << "\t >> m7->Mult(&mout, &m8); // The solution is stored in mout" << endl;
	double datos7[4] = { 1, 2, 3, 1 };
	Matrix *m7 = new Matrix(2, 2, datos7);
	double datos8[4] = { 1, 2, 3, 4 };
	Matrix m8(2, 2, datos8);
	m7->Mult(&mout, &m8);  // The solution is stored in mout
	cout << *m7 << "\n * \n" << m8 << "\n = \n" << mout << endl;

	cout << "\t >> m7->Mult(&m8); // The solution is stored in m7" << endl;
	m7->Mult(&m8);			// The solution is stored in m7
	cout << *m7 << endl;

	cout << "\n*** Matrix ADD: " << endl;
	cout << "\t >> m7->Add(&mout, &m8); // The solution is stored in mout" << endl;
	m7->Add(&mout, &m8);  // The solution is stored in mout
	cout << *m7 << "\n + \n" << m8 << "\n = \n" << mout << endl;

	cout << "\t >> m7->Add(&m8); // The solution is stored in m7" << endl;
	m7->Add(&m8);			// The solution is stored in m7
	cout << *m7 << endl;

	cout << "\n*** Matrix SUBSTRACTION: \n" << endl;
	cout << "\t >> m7->Subs(&mout, &m8); // The solution is stored in mout" << endl;
	m7->Subs(&mout, &m8);  // The solution is stored in mout
	cout << *m7 << "\n - \n" << m8 << "\n = \n" << mout << endl;

	cout << "\t >> m7->Subs(&m8); // The solution is stored in m7" << endl;
	m7->Subs(&m8);			// The solution is stored in m7
	cout << *m7 << endl;

	system("pause");

	cout << "\n*** Matrix MULTIPLICATION by an INTEGER: \n" << endl;
	cout << "\t >> m7->Mult(&mout, 2); // The solution is stored in mout" << endl;
	m7->Mult(&mout, 2);
	cout << *m7 << "\n * \n" << 2 << "\n = \n" << mout << endl;

	cout << "\t >> m7->Mult(2); // The solution is stored in m7" << endl;
	m7->Mult(2);
	cout << *m7 << endl;

	cout << "\n*** Matrix DIVISION by an INTEGER: \n" << endl;
	cout << "\t >> m7->Division(&mout, 2); // The solution is stored in mout" << endl;
	m7->Division(&mout, 2);
	cout << *m7 << "\n / \n" << 2 << "\n = \n" << mout << endl;

	cout << "\t >> m7->Division(2); // The solution is stored in m7" << endl;
	m7->Division(2);
	cout << *m7 << endl;

	system("pause");

	cout << "\n*** Matrix INVERSE3X3: \n" << endl;
	double data9[9] = { 2,3,4,1,2,3,5,6,0 };
	Matrix m9(3, 3, data9);

	cout << "\t >> m9.Inv3x3(&mout); // The solution is stored in mout" << endl;
	m9.Inv3x3(&mout);
	cout << m9 << "\n ^-1 \n" << mout << endl;

	cout << "\t >> m9.Inv3x3(); // The solution is stored in m9" << endl;
	m9.Inv3x3();
	cout << m9 << endl;

	system("pause");

	cout << "\n*** Matrix INVERSE4x4: \n" << endl;
	double dataa9[16] = { 2,3,4,1,2,3,5,6,0,10,1,3,5,6,21,3};
	Matrix mm9(4, 4, dataa9);

	cout << "\t >> mm9.Inv4x4(&mout); // The solution is stored in mout" << endl;
	mm9.Inv4x4(&mout);
	cout << mm9 << "\n ^-1 \n" << mout << endl;

	cout << "\t >> mm9.Inv4x4(); // The solution is stored in mm9" << endl;
	mm9.Inv4x4();
	cout << mm9 << endl;

	system("pause");

	cout << "\n*** Matrix INVERSE GENERIC: \n" << endl;
	cout << "\t >> mm9.Inverse(&mout); // The solution is stored in mout" << endl;
	mm9.Inverse(&mout);
	cout << mm9 << "\n ^-1 \n" << mout << endl;

	cout << "\t >> mm9.Inverse(); // The solution is stored in mm9" << endl;
	mm9.Inverse();
	cout << mm9 << endl;

	system("pause");

	cout << "\n*** Trace of a square Matrix: \n" << endl;
	cout << "\t >> double d;" << endl;
	cout << "\t >> d = m.Trace()" << endl;
	Matrix m11(3, 3, data9);
	double d;
	d = m11.Trace();

	cout << d << endl;

	system("pause");

	return 0;
}

int howTo_Vector()
{

	cout << "\n*** Fill vector with Parenthesis:" << endl;
	cout << "\t >> Vector v(2);" << endl;
	cout << "\t >> v(0) = 1;" << endl;

	Vector v(2);
	v(0) = 1;
	v(1) = 2;
	cout << v << endl;

	cout << "\n*** Print vector: \n" << endl;
	cout << "\t >> cout << v << endl; \n" << endl;

	system("pause");

	cout << "\n*** Fill vector with the constructor: \n" << endl;
	cout << "\t >> double datos[4] = { 1, 2, 3, 4 };" << endl;
	cout << "\t >> Vector v2(2, datos); \n" << endl;

	double datos[4] = { 1, 2, 3, 4 };
	Vector v2(4, datos);

	cout << v2 << endl;

	system("pause");

	cout << "\n*** GET/SET properties: \n" << endl;
	cout << "\t >> v.getLength();" << endl;
	cout << v2.getLength() << endl;

	cout << "\t >> v.setLength(5);" << endl;
	v.setLength(5);
	cout << v << endl;


	system("pause");

	cout << "\n*** STANDARD vectors: \n" << endl;
	cout << "\t >> v.Zeros();" << endl;
	v.Zeros();
	cout << v << endl;

	cout << "\t >> v.Ones();" << endl;
	v.Ones();
	cout << v << endl;

	cout << "\t >> v.Rand(); //TODO" << endl;
	//v2.Rand();
	//cout << v2 << endl;

 	system("pause");
	
	cout << "\n*** TRANSPOSE VECTOR: \n" << endl;
	cout << "\t >> v.Transpose(&v2); // The solution is stored in v2" << endl;
	v.Transpose(&v2);
	cout << "[" << v << "]' = " << endl;
	cout << v2 << endl;

	cout << "\t >> v.Transpose(); // The solution is stored in v" << endl;
	v.Transpose();
	cout << v << endl;


	double datos3[4] = { 1, 2, 3, 4 };
	Vector v3(4, datos3);
	double datos4[4] = { 2, 1, 4, 3 };
	Vector v4(4, datos4);
	Vector outVec;

	cout << "\n*** ADDING VECTOR: \n" << endl;
	cout << "\t >> v3.Add(&outVec, &v4); // The solution is stored in outVec" << endl;
	cout << v3 << endl;
	cout << "+" << endl;
	cout << v4 << endl;
	v3.Add(&outVec, &v4);
	cout << outVec << endl;
	cout << "\t >> v3.Add(&v4); // The solution is stored in v3" << endl;
	v3.Add(&v4);
	cout << v3 << endl;

	cout << "\n*** SUBSTRACTING VECTOR: \n" << endl;
	v3.Transpose();
	v4.Transpose();
	cout << "\t >> v3.Subs(&outVec, &v4); // The solution is stored in outVec" << endl;
	cout << v3 << endl;
	cout << "-" << endl;
	cout << v4 << endl;
	v3.Subs(&outVec, &v4);
	cout << outVec << endl;
	cout << "\t >> v3.Subs(&v4); // The solution is stored in v3" << endl;
	v3.Subs(&v4);
	cout << v3 << endl;

	system("pause");

	cout << "\n*** VECTOR MULTIPLICATION: \n" << endl;
	v3.Transpose();
	cout << "\t >> Vector vout;" << endl;
	cout << "\t >> v3.Mult(&vout, &v4);" << endl;
	cout << v3 << endl;
	cout << "*" << endl;
	cout << v4 << "=\n" << endl;

	Vector vout;
	v3.Mult(&vout, &v4);
	cout << vout << endl;

	if (v3.isRow())
		v3.Transpose();
	if (v4.isCol())
		v4.Transpose();
	cout << "\t >> Matrix mout;" << endl;
	cout << "\t >> v3.Mult(&mout, &v4);" << endl;
	cout << v3 << endl;
	cout << "*" << endl;
	cout << v4 << "=\n" << endl;
	Matrix mout;
	v3.Mult(&mout, &v4);
	cout << mout << endl;

	Vector vout2;
	double data9[9] = { 2,3,4,1,2,3,5,6,0 };
	Matrix m9(3, 3, data9);
	double datav[3] = { 1,2,3 };
	Vector vv(3, datav);
	cout << "\t >> vv.Mult(&vout2, &m9);" << endl;
	vv.Mult(&vout2, &m9);
	cout << vv << endl;
	cout << "*" << endl;
	cout << m9 << "=\n" << endl;
	cout << vout2 << endl;

	cout << "\t >> vv.Mult(&m9);" << endl;
	vv.Mult(&m9);
	cout << vv << endl;

	double dataa[9] = { 2,3,4 };
	Matrix mm(1, 3, dataa);
	double datav1[3] = { 1,2,3 };
	Vector vv1(3, datav1);
	vv1.Transpose();
	cout << "\t >> vv.Mult(&mout, &mm);" << endl;
	vv1.Mult(&mout, &mm);
	cout << vv1 << endl;
	cout << "*" << endl;
	cout << mm << "=\n" << endl;
	cout << mout << endl;
	

	system("pause");

	cout << "\n*** VECTOR ELEMENT MULTIPLICATION: \n" << endl;
	if (v3.isCol())
		v3.Transpose();
	if (v4.isCol())
		v4.Transpose();

	cout << "\t >> v3.ElementMult(&outVec, &v4); // The solution is stored in outVec" << endl;
	cout << v3 << endl;
	cout << ".*" << endl;
	cout << v4 << endl;
	v3.ElementMult(&outVec, &v4);
	cout << outVec << endl;
	cout << "\t >> v3.ElementMult(&v4); // The solution is stored in v3" << endl;
	v3.ElementMult(&v4);
	cout << v3 << endl;

	system("pause");

	cout << "\n*** MULTIPLICATION by an integer: \n" << endl;
	cout << "\t >> v3.Mult(&outVec, 2); // The solution is stored in outVec" << endl;
	cout << v3 << endl;
	cout << "*" << endl;
	cout << 2 << "=\n" << endl;
	v3.Mult(&outVec, 2);
	cout << outVec << endl;
	cout << "\t >> v3.ElementMult(2); // The solution is stored in v3" << endl;
	v3.Mult(2);
	cout << v3 << endl;

	cout << "\n*** DIVISION by an integer: \n" << endl;
	v3.Transpose();
	v4.Transpose();
	cout << "\t >> v3.Division(&outVec, 2); // The solution is stored in outVec" << endl;
	cout << v3 << endl;
	cout << "/" << endl;
	cout << 2 << "=\n" << endl;
	v3.Division(&outVec, 2);
	cout << outVec << endl;
	cout << "\t >> v3.Division(2); // The solution is stored in v3" << endl;
	v3.Division(2);
	cout << v3 << endl;
	
	system("pause");

	cout << "\n*** EXPONENTIAL by an integer: \n" << endl;
	v3.Transpose();
	v4.Transpose();
	cout << "\t >> v3.ElementExp(&outVec, 2); // The solution is stored in outVec" << endl;
	cout << v3 << endl;
	cout << ".^" << endl;
	cout << 2 << "=\n" << endl;
	v3.ElementExp(&outVec, 2);
	cout << outVec << endl;
	cout << "\t >> v3.ElementExp(2); // The solution is stored in v3" << endl;
	v3.ElementExp(2);
	cout << v3 << endl;

	system("pause");

	cout << "\n*** Cross product: \n" << endl;

	double data5[3] = { 1, 2, 3 };
	Vector v5(3, data5);
	double data6[3] = { 3, 2, 1 };
	Vector v6(3, data6);
	Vector v_out;

	cout << "\t >> v5.CrossProd(&v_out, &v6); // The solution is stored in v_out" << endl;
	v5.CrossProd(&v_out, &v6);
	cout << v5 << "\n x \n" << v6 << "=\n" << endl;
	cout << v_out << endl;

	cout << "\t >> v5.CrossProd(&v6);; // The solution is stored in v5" << endl;
	v5.CrossProd(&v6);
	cout << v5 << endl;

	system("pause");

	cout << "\n*** Norm of a vector Nx1 or 1xN: \n" << endl;

	Vector v7(3, data6);
	double norm;
	cout << "\t >> norm = v.Norm();" << endl;
	norm = v7.Norm();
	cout << "Norm([ \n" << v7 << "])\n" << endl;
	cout << norm << endl;

	system("pause");

	cout << "\n*** Vector Normalization: \n" << endl;

	Vector v8(3, data6);

	cout << "\t >> m.VecNorm(&vout); // The solution is stored in vout" << endl;
	v8.VecNorm(&v_out);
	cout << "VecNorm([ \n" << v8 << "])\n" << endl;
	cout << v_out << endl;

	cout << "\t >> v.VecNorm(); // The solution is stored in v" << endl;
	v8.VecNorm();
	cout << v8 << endl;

	system("pause");

	cout << "\n*** Transformations between matrix and vector classes \n" << endl;
	cout << "\t >> vvv.vec2mat(&mv);" << endl;
	Vector vvv(3, data6);
	Matrix mv;
	vvv.Transpose();
	vvv.vec2mat(&mv);
	cout << vvv << endl;
	cout << mv << endl;

	cout << "\t >> vvv.mat2vec(&mv);" << endl;
	mv(0, 0) = 10;
	mv.Transpose();
	vvv.mat2vec(&mv);
	cout << mv << endl;
	cout << vvv << endl;

	system("pause");

	return 1;
}

int howto_Quaternion()
{
	cout << "\n*** THE STRUCTURE OF THE QUATERNION IS:" << endl;
	cout << "q(0) -> qs or qw" << endl;
	cout << "q(1) -> qx" << endl;
	cout << "q(2) -> qy" << endl;
	cout << "q(3) -> qz \n" << endl;
	cout << "\n*** Fill Quaternion with Parenthesis:" << endl;
	cout << "\t >> Quaternion q();" << endl;
	cout << "\t >> q(0) = 1;" << endl;
	Quaternion q;
	q(0) = 1;
	q(1) = 2;
	q(2) = 3;
	q(3) = 4;
	cout << q << endl;

	cout << "\n*** Fill matrix with the constructor: " << endl;
	cout << "\t Quaternion::QuaternionStruct q_struct;" << endl;
	cout << "\t >> q_struct.qs = 0.8; ..." << endl;
	cout << "\t >> Quaternion q2(datos);" << endl;
	Quaternion::QuaternionStruct q_struct;
	q_struct.qs = 0.8;
	q_struct.qx = 0.6;
	q_struct.qy = 1;
	q_struct.qz = -0.2;
	Quaternion q2(&q_struct);
	cout << q2 << endl;

	cout << "\t >> Quaternion q3(0.4, 0.8, -0.7, 0.5);" << endl;
	Quaternion q3(0.4, 0.8, -0.7, 0.5);
	cout << q3 << endl;

	cout << "\n*** Print Quaternion:" << endl;
	cout << "\t >> cout << q << endl; \n" << endl;

	system("pause");

	cout << "\n*** Norm Quaternion:" << endl;
	cout << "\t >> q.Norm(&qout) // The solution is stored in qout;\n" << endl;
	Quaternion qout;
	q.Norm(&qout);
	
	cout << qout << endl;

	cout << "\t >> q.Norm() // The solution is stored in q;\n" << endl;
	q.Norm();
	cout << q << endl;

	cout << "\n*** Quaternion Product:" << endl;
	cout << "\t >> q.Prod(&qout,q2) // The solution is stored in qout;\n" << endl;
	q.Prod(&qout, &q2);

	cout << qout << endl;

	cout << "\t >> q.Prod(q2) // The solution is stored in q;\n" << endl;
	q.Prod(&q2);
	cout << q << endl;

	cout << "\n*** Quaternion Conj:" << endl;
	cout << "\t >> q.Conj(&qout) // The solution is stored in qout;\n" << endl;
	q.Conj(&qout);

	cout << qout << endl;

	cout << "\t >> q.Conj() // The solution is stored in q;\n" << endl;
	q.Conj();
	cout << q << endl;

	system("pause");
	cout << "\n*** Quaternion to Homogeneous Matrix:" << endl;
	cout << "\t >> Matrix m;" << endl;
	cout << "\t >> m = q.q2tr();\n" << endl;
	Matrix m;
	m = q.q2tr();

	cout << m << endl;

	cout << "\n*** Quaternion to Rotation Matrix:" << endl;
	cout << "\t >> Matrix m;" << endl;
	cout << "\t >> m = q.q2r();\n" << endl;
	m = q.q2r();

	cout << m << endl;

	cout << "\n*** Quaternion to Euler Angles (Roll-Pitch-Yaw)(XYZ):" << endl;
	cout << "\t >> Vector v;" << endl;
	cout << "\t >> v = q.q2Euler();\n" << endl;
	Vector v;
	Quaternion q_euler;
	q(0) = -0.9948;
	q(1) =  0.0124;
	q(2) = 0.0042;
	q(3) = 0.1013;
	v = q.q2Euler();

	cout << v << endl;

	system("pause");

	cout << "\n*** Quaternion Update:" << endl;
	cout << "\t Quaternion::QuaternionStruct q_struct;" << endl;
	cout << "\t >> q_struct.qs = 0.8; ..." << endl;
	cout << "\t >> q.Update(&q_struct);" << endl;
	Quaternion::QuaternionStruct q_struct2;
	q_struct2.qs = 0.8;
	q_struct2.qx = 0.6;
	q_struct2.qy = 1;
	q_struct2.qz = -0.2;
	
	q.Update(&q_struct2);
	cout << q << endl;
	cout << "\t >> q.Update(0.4,0.8,-0.7,0.5);" << endl;
	q.Update(0.4,0.8,-0.7,0.5);
	cout << q << endl;
	
	system("pause");

	return 0;
}

int main()
{
	cout.setf(std::ios::fixed, std::ios::floatfield);
	cout.precision(3);

	cout << "HOW TO USE MATRIX CLASS OF THE LIBRARY SIMPLE MATH" << endl;
	howTo_Matrix();
	cout << "END OF MATRIX CLASS" << endl;
	system("pause");

	system("CLS");

	cout << "HOW TO USE VECTOR CLASS OF THE LIBRARY SIMPLE MATH" << endl;
	howTo_Vector();
	cout << "END OF VECTOR CLASS" << endl;
	system("pause");

	system("CLS");
	
	cout << "HOW TO USE QUATERNION CLASS OF THE LIBRARY SIMPLE MATH" << endl;
	howto_Quaternion();
	cout << "END OF QUATERNION CLASS" << endl;
	system("pause");

	return 1;
}

