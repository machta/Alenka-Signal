#ifndef WX_H
#define	WX_H

#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <memory>
#include <cassert>

template<class T>
using wxVector = std::vector<T>;

#define wxT(a_) wxString(a_)

const int wxEVT_THREAD = 0;

class wxString : public std::string
{
public:
	wxString() : std::string() {}
	wxString(const char* str) : std::string(str) {}

	template<class... T>
	int Printf(const wxString& pszFormat, T... p)
	{
		const char* format = pszFormat.c_str();

		int size1 = snprintf(nullptr, 0, format, p...);

		std::unique_ptr<char[]> buffer(new char[size1 + 1]);

		int size2 = sprintf(buffer.get(), format, p...);
		assert(size1 == size2); (void)size2;

		assign(buffer.get(), buffer.get() + size1);

		return size1;
	}
};

class wxThreadEvent
{
	int intVal;
	wxString stringVal;

public:
	wxThreadEvent(int eventType, int id) { (void)eventType; (void)id;}

	void SetInt(int intCommand)
	{
		intVal = intCommand;
	}

	void SetString(const wxString& string)
	{
		stringVal = string;
	}

	wxThreadEvent* Clone() const
	{
		return new wxThreadEvent(0, 0);
	}
};

class wxEvtHandler
{};

inline void wxQueueEvent(wxEvtHandler* dest, wxThreadEvent* event)
{
	(void)dest; (void)event;
}

#endif // WX_H
