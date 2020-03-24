#ifndef AtcRegion_h
# define AtcRegion_h

# include <vector>

class AtcBox {         //������� ��������
public:
// ���������� ����� ��������: z - � ���������, r - � ��������
	float lMin, lMax;
// ������������ ���� ����� �� ��������� �������� (������ - ����������)
	float phiMin, phiMax;
	AtcBox(float l1, float phi1, float l2, float phi2) :       //������� ��������
		lMin(l1<?l2),
		lMax(l1>?l2),
		phiMin(phi1<?phi2),
		phiMax(phi1>?phi2)
	{}
};

class AtcRegion     //������� ��������
{
public:
	std::vector<AtcBox> region;

	void addBox(float l1, float phi1, float l2, float phi2) {     //������� ��������� ������ ����� ������� ��������
		region.push_back(AtcBox(l1,phi1,l2,phi2));
	}
	void addBoxD(float l1, float phi1, float l2, float phi2) {     //������� ��� ������� �������� (���������)
		region.push_back(AtcBox(l1,phi1,l2,phi2));
		region.push_back(AtcBox(l1,-phi2,l2,-phi1));
	}

	bool isIn(float x, float phi) const;
	bool isIn(float x1, float phi1, float x2, float phi2) const;
};

inline bool AtcRegion::isIn(float x, float phi) const
{
	std::vector<AtcBox>::const_iterator iter=region.begin(), last=region.end();

	for(; iter!=last; iter++) {
		const AtcBox& box=*iter;
		if( x>=box.lMin && x<=box.lMax && phi>=box.phiMin && phi<=box.phiMax )
			return true;
	}
	return false;
}
inline bool AtcRegion::isIn(float x1, float phi1, float x2, float phi2) const
{
	return (isIn(x1,phi1) && isIn(x2,phi2));
}

#endif //AtcRegion_h
