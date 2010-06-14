#pragma once
namespace Another
{
	template<class Value>
	class Sort
	{
	public:
		bool operator () ( Value const &_A, Value const &_B) const
		{
			if(_A->GetIndex() <= _B->GetIndex() ) return true;
			return false;
		}
	};

}
