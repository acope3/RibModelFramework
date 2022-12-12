
// This is a mapping of delta category to mixture definition.
class mixtureDefinition
{
public:
	unsigned delM; // With PA(NSE) model, this refers to alpha
	unsigned delEta; //With PA(NSE) model, this refers to lambda
	unsigned nse; //With PANSE, this refers to nse
	unsigned phi; //Previous versions have linked selection and synthesis categories, but we would like the option to decouple this in the case of PA/PANSE

};
