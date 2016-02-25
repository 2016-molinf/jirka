$(document).ready(function()
{
    function getCookie(name)
    {
    var cookieValue = null;
    if (document.cookie && document.cookie != '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = jQuery.trim(cookies[i]);
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) == (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
    }

    var csrftoken = getCookie('csrftoken');

    function csrfSafeMethod(method)
    {
        return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
    }

    $.ajaxSetup({
        beforeSend: function(xhr, settings)
        {
            if (!csrfSafeMethod(settings.type) && !this.crossDomain)
            {
                xhr.setRequestHeader("X-CSRFToken", csrftoken);
            }
        }
});

    var mol = ChemDoodle.readMOL(`C2H6
APtclcactv02241615052D 0   0.00000     0.00000

  8  7  0  0  0  0  0  0  0  0999 V2000
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.3800    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -0.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000   -0.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.6200   -0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
M  END
$$$$`);

    console.log(mol);

    $("#button-test").click(function() {
        $.ajax({
            url: "/api/api_molConverter",
            method: "POST",
            data: {data: "CC", format_from: "smiles", format_to: "molfile"},
            success: function(result) {
                if (result.success)
                {
                    var mol = ChemDoodle.readMOL(result.data);
                    console.log(mol);
                    console.log(document.sketcher);
                    sketcher.loadMolecule(mol);
                }
            }
        });
    });
});