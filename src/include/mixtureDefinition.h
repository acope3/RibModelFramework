
// This is a mapping of delta category to mixture definition.
class mixtureDefinition
{
public:
	int delM = -1; // With PA(NSE) model, this refers to alpha
	int delEta = -1; //With PA(NSE) model, this refers to lambda
	int nse = -1; //With PANSE, this refers to nse
	int phi = -1; //Previous versions have linked selection and synthesis categories, but we would like the option to decouple this in the case of PA/PANSE

};
