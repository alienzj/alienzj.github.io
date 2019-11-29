let image = document.querySelector('img');

image.onclick = function() {
    let src = image.getAttribute('src');
    if (src == 'images/hello_world.jpg') {
        image.setAttribute('src', 'images/dna.jpg');
    } else {
        image.setAttribute('src', 'images/hello_world.jpg');
    }
}

let button = document.querySelector('button');
let heading = document.querySelector('h1');

function set_user_name() {
    let name = prompt("Please input your name: ");
    if (!name || name === null) {
        set_user_name();
    } else {
        localStorage.setItem('name', name);
        heading.textContent = 'Hello, ' + name;
    }
}

if (!localStorage.getItem('name')) {
    set_user_name();
} else {
    let store_name = localStorage.getItem('name');
    heading.textContent = 'Hello, ' + store_name;
}

button.onclick = function() {
    set_user_name();
}
