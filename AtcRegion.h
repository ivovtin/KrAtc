#ifndef AtcRegion_h
# define AtcRegion_h

# include <vector>

class AtcBox {         //коробка счетчика
public:
// координаты вдоль счетчика: z - в барельном, r - в торцевом
	float lMin, lMax;
// азимутальные углы одной из половинок счетчика (вторая - зеркальная)
	float phiMin, phiMax;
	AtcBox(float l1, float phi1, float l2, float phi2) :       //размеры счетчика
		lMin(l1<?l2),
		lMax(l1>?l2),
		phiMin(phi1<?phi2),
		phiMax(phi1>?phi2)
	{}
};

class AtcRegion     //область счетчика
{
public:
	std::vector<AtcBox> region;

	void addBox(float l1, float phi1, float l2, float phi2) {     //смотрим половинку только одной области счетчика
		region.push_back(AtcBox(l1,phi1,l2,phi2));
	}
	void addBoxD(float l1, float phi1, float l2, float phi2) {     //смотрим две области счетчика (зеркально)
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
