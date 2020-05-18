Vue.component("menu-comp", {
	template: `
	<div @mouseover="active=true" @mouseleave="active=false">
		<div id="contents-button">
			<button>Menu</button>
		</div>
	    <div id="contents" v-if="active">
	    	<ul id="content-list"> 
	         	<li v-for="element in elements"> <a :href="element.htmlRef" class="content-links"> {{ element.name }} </a></li>
	        </ul> 
	    </div>
	</div>
	`,
	data() {
		return {
			elements: [
				{
					id: 1,
					name: "Home",
					htmlRef: "index.html"
				},
				{
					id: 2,
					name: "Logg",
					htmlRef: "logg.html"
				},
				{
					id: 3,
					name: "JavaScript",
					htmlRef: "pagejavascript.html"
				}],
			active: false		
		}
	},
	methods: {
		stuff() {
			console.log("Stuff")
			this.active = !this.active
		}
	}
})

Vue.component("post-input", {
	template: `
	<form class="blogPost" @submit.prevent="onSubmit">
		<p>
			<label><h2>New post:</h2></label>
			<textarea></textarea>
			<br>
			<input type="submit" value="Submit">
		</p>
	</form>
	`,
	data() {
		return {
			posts: []
		}
	},
	methods: {
		onSubmit() {
			this.posts.push("post-added")
		}
	}
})

var page = new Vue({
	el: "#page",
})

